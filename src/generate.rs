//! Reproducible random sequence generator for benchmarking.
//!
//! Ambiguous bases can be placed two ways:
//!   * **scattered** (`*_runs == 0`): each position independently becomes the
//!     code with the given probability — a stress test that yields millions of
//!     length-1 blocks.
//!   * **clustered** (`*_runs == K`): the budget (`frac * length`) is split into
//!     K random-sized runs placed at random disjoint positions — like the long
//!     assembly-gap `N` runs (and ambiguous regions) of a real genome.
//!
//! Uses SplitMix64 so test data is identical across runs/platforms for a seed.

use crate::seq::Sequence;

/// SplitMix64 PRNG (public-domain algorithm by Sebastiano Vigna).
pub struct SplitMix64(u64);

impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        SplitMix64(seed)
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform float in [0, 1).
    pub fn frac(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform integer in [0, n).
    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

pub struct GenOpts {
    pub seqs: usize,
    pub length: usize,
    pub n_frac: f64,
    pub iub_frac: f64,
    pub seed: u64,
    pub name_prefix: String,
    /// 0 => scatter N per-base; K => K clustered N runs per sequence.
    pub n_runs: usize,
    /// 0 => scatter IUB per-base; K => K clustered IUB runs per sequence.
    pub iub_runs: usize,
}

const ACGT: &[u8; 4] = b"ACGT";
/// The ten specific IUB degenerate codes (upper-case only, per the test spec).
const IUB: &[u8; 10] = b"RYSWKMBDHV";

/// Split `total` into `parts` positive integers summing to `total` (random).
fn random_composition(rng: &mut SplitMix64, total: usize, parts: usize) -> Vec<usize> {
    if parts == 0 || total == 0 {
        return Vec::new();
    }
    let parts = parts.min(total);
    let mut sizes = vec![1usize; parts];
    let mut remaining = total - parts;
    for s in sizes.iter_mut().take(parts - 1) {
        let take = rng.below(remaining + 1);
        *s += take;
        remaining -= take;
    }
    sizes[parts - 1] += remaining;
    sizes
}

/// Split `free` into `slots` non-negative integers summing to `free` (random).
fn random_gaps(rng: &mut SplitMix64, free: usize, slots: usize) -> Vec<usize> {
    let mut gaps = vec![0usize; slots];
    if slots == 0 {
        return gaps;
    }
    let mut remaining = free;
    for g in gaps.iter_mut().take(slots - 1) {
        let take = rng.below(remaining + 1);
        *g = take;
        remaining -= take;
    }
    gaps[slots - 1] = remaining;
    gaps
}

/// Generate `opts.seqs` sequences with the requested ambiguity content.
pub fn generate(opts: &GenOpts) -> Vec<Sequence> {
    let mut rng = SplitMix64::new(opts.seed);
    let mut out = Vec::with_capacity(opts.seqs);
    let len = opts.length;

    for s in 0..opts.seqs {
        let mut bases: Vec<u8> = (0..len).map(|_| ACGT[rng.below(4)]).collect();

        // Scattered ambiguity (only for kinds not being clustered).
        if opts.iub_runs == 0 {
            for b in bases.iter_mut() {
                if rng.frac() < opts.iub_frac {
                    *b = IUB[rng.below(IUB.len())];
                }
            }
        }
        if opts.n_runs == 0 {
            for b in bases.iter_mut() {
                if rng.frac() < opts.n_frac {
                    *b = b'N';
                }
            }
        }

        // Clustered ambiguity: build run sizes for each kind, then lay all runs
        // down at random disjoint positions. kind 0 = N, kind 1 = IUB.
        let mut runs: Vec<(usize, u8)> = Vec::new();
        if opts.n_runs > 0 {
            let total = (opts.n_frac * len as f64).round() as usize;
            runs.extend(random_composition(&mut rng, total, opts.n_runs).into_iter().map(|sz| (sz, 0)));
        }
        if opts.iub_runs > 0 {
            let total = (opts.iub_frac * len as f64).round() as usize;
            runs.extend(random_composition(&mut rng, total, opts.iub_runs).into_iter().map(|sz| (sz, 1)));
        }
        if !runs.is_empty() {
            // Fisher–Yates shuffle so N and IUB runs interleave randomly.
            for i in (1..runs.len()).rev() {
                let j = rng.below(i + 1);
                runs.swap(i, j);
            }
            let used: usize = runs.iter().map(|(sz, _)| *sz).sum();
            let gaps = random_gaps(&mut rng, len.saturating_sub(used), runs.len() + 1);
            let mut pos = 0usize;
            for (k, &(sz, kind)) in runs.iter().enumerate() {
                pos += gaps[k];
                let end = (pos + sz).min(len);
                for j in pos..end {
                    bases[j] = if kind == 0 { b'N' } else { IUB[rng.below(IUB.len())] };
                }
                pos = end;
            }
        }

        out.push(Sequence::new(format!("{}{}", opts.name_prefix, s), bases));
    }
    out
}
