//! A minimal on-disk B+ tree mapping a sequence name → u64 value (a record
//! offset), modelled on the UCSC bigBed `bPlusTree`: fixed-width keys, fan-out
//! `BLOCK_SIZE`, O(log N) lookup with no full load.
//!
//! The whole tree is serialized into one position-independent blob: every node
//! reference is an offset **relative to the start of the blob**, so the blob can
//! be dropped anywhere in a file and read with `Bpt::new(&file[tocOffset..])`.
//! Leaf *values*, by contrast, are stored verbatim (we use absolute file
//! offsets to the sequence records).
//!
//! ```text
//! blob header (24 bytes, little-endian)
//!   keySize    u32        fixed key width = longest name
//!   blockSize  u32        max children/items per node
//!   itemCount  u64
//!   rootOffset u64        blob-relative offset of the root node
//! node
//!   isLeaf     u8 (1/0)
//!   pad        u8
//!   count      u16
//!   leaf:     count × ( key[keySize], value:u64 )
//!   internal: count × ( key[keySize], childOffset:u64 )   key = subtree min key
//! ```

use crate::io::{peek_u16, peek_u32, peek_u64, put_u16, put_u32, put_u64};

const BLOCK_SIZE: usize = 256;
const HDR: usize = 24;

enum Node {
    Leaf(Vec<(Vec<u8>, u64)>),
    Internal(Vec<(Vec<u8>, usize)>), // (subtree-min key, child node index)
}

fn first_key(n: &Node) -> &[u8] {
    match n {
        Node::Leaf(v) => &v[0].0,
        Node::Internal(v) => &v[0].0,
    }
}

/// Serialize a name→value map into a B+ tree blob.
pub fn build(mut items: Vec<(String, u64)>) -> Vec<u8> {
    items.sort_by(|a, b| a.0.as_bytes().cmp(b.0.as_bytes()));
    let key_size = items.iter().map(|(n, _)| n.len()).max().unwrap_or(1).max(1);
    let pad = |s: &str| {
        let mut k = s.as_bytes().to_vec();
        k.resize(key_size, 0);
        k
    };

    // Build bottom-up: a level of leaves, then internal levels until one root.
    let mut nodes: Vec<Node> = Vec::new();
    let mut level: Vec<usize> = Vec::new();
    if items.is_empty() {
        nodes.push(Node::Leaf(Vec::new()));
        level.push(0);
    } else {
        for chunk in items.chunks(BLOCK_SIZE) {
            nodes.push(Node::Leaf(chunk.iter().map(|(n, v)| (pad(n), *v)).collect()));
            level.push(nodes.len() - 1);
        }
    }
    while level.len() > 1 {
        let mut next = Vec::new();
        for chunk in level.chunks(BLOCK_SIZE) {
            let children = chunk.iter().map(|&ci| (first_key(&nodes[ci]).to_vec(), ci)).collect();
            nodes.push(Node::Internal(children));
            next.push(nodes.len() - 1);
        }
        level = next;
    }
    let root = level[0];

    // Assign blob-relative offsets, then serialize.
    let node_size = |n: &Node| -> usize {
        let count = match n {
            Node::Leaf(v) => v.len(),
            Node::Internal(v) => v.len(),
        };
        4 + count * (key_size + 8)
    };
    let mut offsets = vec![0u64; nodes.len()];
    let mut cur = HDR as u64;
    for (i, n) in nodes.iter().enumerate() {
        offsets[i] = cur;
        cur += node_size(n) as u64;
    }

    let mut out = Vec::with_capacity(cur as usize);
    put_u32(&mut out, key_size as u32);
    put_u32(&mut out, BLOCK_SIZE as u32);
    put_u64(&mut out, items.len() as u64);
    put_u64(&mut out, offsets[root]);
    for n in &nodes {
        match n {
            Node::Leaf(v) => {
                out.push(1);
                out.push(0);
                put_u16(&mut out, v.len() as u16);
                for (k, val) in v {
                    out.extend_from_slice(k);
                    put_u64(&mut out, *val);
                }
            }
            Node::Internal(v) => {
                out.push(0);
                out.push(0);
                put_u16(&mut out, v.len() as u16);
                for (k, ci) in v {
                    out.extend_from_slice(k);
                    put_u64(&mut out, offsets[*ci]);
                }
            }
        }
    }
    out
}

/// Reader over a B+ tree blob.
pub struct Bpt<'a> {
    blob: &'a [u8],
    key_size: usize,
    root: usize,
}

impl<'a> Bpt<'a> {
    pub fn new(blob: &'a [u8]) -> Option<Self> {
        if blob.len() < HDR {
            return None;
        }
        let key_size = peek_u32(blob, 0, true)? as usize;
        let root = peek_u64(blob, 16, true)? as usize;
        Some(Bpt { blob, key_size, root })
    }

    pub fn item_count(&self) -> u64 {
        peek_u64(self.blob, 8, true).unwrap_or(0)
    }

    fn g16(&self, o: usize) -> usize {
        peek_u16(self.blob, o, true).unwrap_or(0) as usize
    }
    fn g64(&self, o: usize) -> u64 {
        peek_u64(self.blob, o, true).unwrap_or(0)
    }

    /// Look up a name. O(log N), touching only the nodes on the path.
    pub fn find(&self, name: &str) -> Option<u64> {
        if name.len() > self.key_size {
            return None;
        }
        let mut q = name.as_bytes().to_vec();
        q.resize(self.key_size, 0);
        let stride = self.key_size + 8;
        let mut off = self.root;
        loop {
            let is_leaf = self.blob[off] == 1;
            let count = self.g16(off + 2);
            let entries = off + 4;
            if is_leaf {
                for i in 0..count {
                    let base = entries + i * stride;
                    if &self.blob[base..base + self.key_size] == q.as_slice() {
                        return Some(self.g64(base + self.key_size));
                    }
                }
                return None;
            }
            // Internal: descend into the rightmost child whose min-key <= q.
            let mut child = self.g64(entries + self.key_size); // child 0 by default
            for i in 0..count {
                let base = entries + i * stride;
                if &self.blob[base..base + self.key_size] <= q.as_slice() {
                    child = self.g64(base + self.key_size);
                } else {
                    break;
                }
            }
            off = child as usize;
        }
    }

    /// All (name, value) pairs, in sorted key order (for listing / iteration).
    pub fn iter_all(&self) -> Vec<(String, u64)> {
        let mut out = Vec::new();
        self.walk(self.root, &mut out);
        out
    }

    fn walk(&self, off: usize, out: &mut Vec<(String, u64)>) {
        let is_leaf = self.blob[off] == 1;
        let count = self.g16(off + 2);
        let entries = off + 4;
        let stride = self.key_size + 8;
        for i in 0..count {
            let base = entries + i * stride;
            let val = self.g64(base + self.key_size);
            if is_leaf {
                let key = &self.blob[base..base + self.key_size];
                let end = key.iter().position(|&b| b == 0).unwrap_or(self.key_size);
                out.push((String::from_utf8_lossy(&key[..end]).into_owned(), val));
            } else {
                self.walk(val as usize, out);
            }
        }
    }
}
