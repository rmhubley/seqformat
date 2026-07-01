//! `seqformat` command-line tool.

use seqformat::error::{Error, Result};
use seqformat::generate::{generate, GenOpts};
use seqformat::{fasta, fourbit, samtools, twobit, twobyte};
use std::collections::{HashMap, HashSet};
use std::process::ExitCode;

const USAGE: &str = "\
seqformat — twoBit / 4-bit sequence tool

USAGE:
    seqformat <command> [options]

COMMANDS:
    fa2twobit  <in.fa> <out.2bit> [--long] [--iub] [--index]
        Convert FASTA to twoBit.
          --long   write 64-bit index offsets (UCSC '-long' format)
          --iub    record exact IUB degenerate codes (backward-compatible)
          --index  append a backward-compatible sorted-name index for O(log N)
                   lookup (a footer of 8 bytes/seq; old readers ignore it)

    twobit2fa  <in.2bit> <out.fa> [--width N]
        Convert twoBit (standard / long / IUB-extended) to FASTA.

    fa2fourbit <in.fa> <out.4bit>
    fourbit2fa <in.4bit> <out.fa> [--width N]
        FASTA <-> BWA/BAM-style 4-bit (case not preserved).

    fa2faidx   <in.fa> <out.fa> [--bgzip] [--width N] [--level N]
        Write a samtools-style indexed FASTA: out.fa + out.fa.fai, and with
        --bgzip a BGZF-compressed out.fa.gz + .fai + .gzi (level default 6).

    fa2be      <in.fa> <out.2be>
    be2fa      <in.2be> <out.fa> [--width N]
        Convert FASTA to/from the experimental 2be format: B+ tree TOC +
        per-sequence merged tagged-edit stream (N runs, IUB points/runs, mask).

    extract    <file.2bit|file.4bit> [region | --seq-list <file>] [--out <fa>] [--width N]
        Random-access extraction. A region is  name  or  name:start-end
        (0-based, half-open — same syntax as twoBitToFa -seqList).
        Without a region or --seq-list, extracts every sequence.

    random     <out.fa> [--seqs K] [--length L] [--n-frac F] [--iub-frac F]
                        [--seed S] [--width N] [--prefix P]
                        [--n-runs K] [--iub-runs K]
        Generate reproducible random uppercase test data (default:
        --seqs 1 --length 1000000 --n-frac 0.01 --iub-frac 0.005 --seed 1).
        --n-runs/--iub-runs K cluster that ambiguity into K random-sized runs
        per sequence (assembly-gap style) instead of scattering it (0 = scatter).

    info       <file>            Print format + per-sequence statistics.
    help                         Show this message.
";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(args: &[String]) -> Result<()> {
    let cmd = match args.first() {
        Some(c) => c.as_str(),
        None => {
            print!("{USAGE}");
            return Ok(());
        }
    };
    let rest = &args[1..];
    match cmd {
        "help" | "-h" | "--help" => {
            print!("{USAGE}");
            Ok(())
        }
        "fa2twobit" => cmd_fa2twobit(rest),
        "twobit2fa" => cmd_twobit2fa(rest),
        "fa2fourbit" => cmd_fa2fourbit(rest),
        "fourbit2fa" => cmd_fourbit2fa(rest),
        "fa2faidx" => cmd_fa2faidx(rest),
        "fa2be" => cmd_fa2be(rest),
        "be2fa" => cmd_be2fa(rest),
        "extract" => cmd_extract(rest),
        "random" => cmd_random(rest),
        "info" => cmd_info(rest),
        other => fail(format!("unknown command '{other}' (try 'seqformat help')")),
    }
}

fn fail<T>(msg: impl Into<String>) -> Result<T> {
    Err(Error::Format(msg.into()))
}

// ---- generic option parsing ----

struct Opts {
    pos: Vec<String>,
    flags: HashSet<String>,
    vals: HashMap<String, String>,
}

impl Opts {
    fn val_usize(&self, key: &str, default: usize) -> Result<usize> {
        match self.vals.get(key) {
            Some(v) => v
                .parse()
                .map_err(|_| Error::Format(format!("invalid value for --{key}: '{v}'"))),
            None => Ok(default),
        }
    }
    fn val_f64(&self, key: &str, default: f64) -> Result<f64> {
        match self.vals.get(key) {
            Some(v) => v
                .parse()
                .map_err(|_| Error::Format(format!("invalid value for --{key}: '{v}'"))),
            None => Ok(default),
        }
    }
    fn val_u64(&self, key: &str, default: u64) -> Result<u64> {
        match self.vals.get(key) {
            Some(v) => v
                .parse()
                .map_err(|_| Error::Format(format!("invalid value for --{key}: '{v}'"))),
            None => Ok(default),
        }
    }
}

/// `bool_flags` are valueless; `val_flags` consume the next argument.
fn parse(args: &[String], bool_flags: &[&str], val_flags: &[&str]) -> Result<Opts> {
    let mut o = Opts {
        pos: Vec::new(),
        flags: HashSet::new(),
        vals: HashMap::new(),
    };
    let mut it = args.iter();
    while let Some(a) = it.next() {
        if let Some(name) = a.strip_prefix("--") {
            if bool_flags.contains(&name) {
                o.flags.insert(name.to_string());
            } else if val_flags.contains(&name) {
                let v = it
                    .next()
                    .ok_or_else(|| Error::Format(format!("--{name} needs a value")))?;
                o.vals.insert(name.to_string(), v.clone());
            } else {
                return fail(format!("unknown option '--{name}'"));
            }
        } else {
            o.pos.push(a.clone());
        }
    }
    Ok(o)
}

// ---- commands ----

fn cmd_fa2twobit(args: &[String]) -> Result<()> {
    let o = parse(args, &["long", "iub", "index"], &[])?;
    let (input, output) = two(&o, "fa2twobit")?;
    let seqs = fasta::read_file(input)?;
    let long = o.flags.contains("long");
    let iub = o.flags.contains("iub");
    let index = o.flags.contains("index");
    if index {
        twobit::write_file_indexed(output, &seqs, long, iub)?;
    } else {
        twobit::write_file(output, &seqs, long, iub)?;
    }
    eprintln!(
        "wrote {} sequence(s) to {output} (twoBit {}{}{})",
        seqs.len(),
        if long { "long" } else { "standard" },
        if iub { " + IUB extension" } else { "" },
        if index { " + name index" } else { "" }
    );
    Ok(())
}

fn cmd_twobit2fa(args: &[String]) -> Result<()> {
    let o = parse(args, &[], &["width"])?;
    let (input, output) = two(&o, "twobit2fa")?;
    let tb = twobit::read_file(input)?;
    fasta::write_file(output, &tb.sequences, o.val_usize("width", 60)?)?;
    eprintln!(
        "wrote {} sequence(s) to {output} (twoBit {}{}{})",
        tb.sequences.len(),
        if tb.long { "long" } else { "standard" },
        if tb.iub { " + IUB extension" } else { "" },
        if tb.indexed { " + name index" } else { "" }
    );
    Ok(())
}

fn cmd_fa2fourbit(args: &[String]) -> Result<()> {
    let o = parse(args, &[], &[])?;
    let (input, output) = two(&o, "fa2fourbit")?;
    let seqs = fasta::read_file(input)?;
    fourbit::write_file(output, &seqs)?;
    eprintln!("wrote {} sequence(s) to {output} (4-bit)", seqs.len());
    Ok(())
}

fn cmd_fourbit2fa(args: &[String]) -> Result<()> {
    let o = parse(args, &[], &["width"])?;
    let (input, output) = two(&o, "fourbit2fa")?;
    let seqs = fourbit::read_file(input)?;
    fasta::write_file(output, &seqs, o.val_usize("width", 60)?)?;
    eprintln!("wrote {} sequence(s) to {output} (read 4-bit)", seqs.len());
    Ok(())
}

fn cmd_fa2faidx(args: &[String]) -> Result<()> {
    let o = parse(args, &["bgzip"], &["width", "level"])?;
    let (input, output) = two(&o, "fa2faidx")?;
    let seqs = fasta::read_file(input)?;
    let bgzip = o.flags.contains("bgzip");
    samtools::write_file(
        output,
        &seqs,
        o.val_usize("width", 60)?,
        bgzip,
        o.val_u64("level", 6)? as u32,
    )?;
    eprintln!(
        "wrote {} sequence(s) to {output}{} (samtools indexed FASTA{})",
        seqs.len(),
        if bgzip { ".fai/.gzi" } else { " + .fai" },
        if bgzip { ", BGZF" } else { "" }
    );
    Ok(())
}

fn cmd_fa2be(args: &[String]) -> Result<()> {
    let o = parse(args, &[], &[])?;
    let (input, output) = two(&o, "fa2be")?;
    let seqs = fasta::read_file(input)?;
    twobyte::write_file(output, &seqs)?;
    eprintln!("wrote {} sequence(s) to {output} (2be)", seqs.len());
    Ok(())
}

fn cmd_be2fa(args: &[String]) -> Result<()> {
    let o = parse(args, &[], &["width"])?;
    let (input, output) = two(&o, "be2fa")?;
    let seqs = twobyte::TwoByteReader::open(input)?.read_all()?;
    fasta::write_file(output, &seqs, o.val_usize("width", 60)?)?;
    eprintln!("wrote {} sequence(s) to {output} (read 2be)", seqs.len());
    Ok(())
}

/// Parse a region spec `name` or `name:start-end` (0-based half-open).
fn parse_region(spec: &str) -> Result<(String, usize, Option<usize>)> {
    match spec.rsplit_once(':') {
        Some((name, range)) => {
            let (s, e) = range
                .split_once('-')
                .ok_or_else(|| Error::Format(format!("bad region '{spec}' (want name:start-end)")))?;
            let start: usize = s
                .parse()
                .map_err(|_| Error::Format(format!("bad start in '{spec}'")))?;
            let end: usize = e
                .parse()
                .map_err(|_| Error::Format(format!("bad end in '{spec}'")))?;
            Ok((name.to_string(), start, Some(end)))
        }
        None => Ok((spec.to_string(), 0, None)),
    }
}

fn region_name(name: &str, start: usize, end: Option<usize>) -> String {
    match end {
        Some(e) => format!("{name}:{start}-{e}"),
        None => name.to_string(),
    }
}

fn cmd_extract(args: &[String]) -> Result<()> {
    let o = parse(args, &["http-stats"], &["out", "width", "seq-list"])?;
    let file = o
        .pos
        .first()
        .ok_or_else(|| Error::Format("extract needs an input file".into()))?;
    let width = o.val_usize("width", 60)?;

    // Build the list of regions to extract.
    let mut regions: Vec<(String, usize, Option<usize>)> = Vec::new();
    if let Some(list) = o.vals.get("seq-list") {
        for line in std::fs::read_to_string(list)?.lines() {
            let line = line.trim();
            if !line.is_empty() {
                regions.push(parse_region(line)?);
            }
        }
    }
    for spec in o.pos.iter().skip(1) {
        regions.push(parse_region(spec)?);
    }

    // Dispatch on the detected format. Each reader exposes the same
    // names()/extract(name,start,end) shape, so a small macro keeps it DRY.
    macro_rules! run {
        ($rd:expr) => {{
            let rd = $rd;
            if regions.is_empty() {
                regions = rd.names().iter().map(|n| (n.clone(), 0, None)).collect();
            }
            regions
                .iter()
                .map(|(n, s, e)| {
                    rd.extract(n, *s, *e)
                        .map(|b| seqformat::Sequence::new(region_name(n, *s, *e), b))
                })
                .collect::<Result<_>>()?
        }};
    }

    // A remote http(s) input is opened over HTTP range requests (UDC-style),
    // never slurped. Format is auto-detected from a small prefix. This is the
    // web-serving path the webseq benchmark exercises.
    if file.starts_with("http://") || file.starts_with("https://") {
        let rd = seqformat::open_url(file)?;
        if regions.is_empty() {
            regions = rd.names().iter().map(|n| (n.clone(), 0, None)).collect();
        }
        let out_seqs: Vec<seqformat::Sequence> = regions
            .iter()
            .map(|(n, s, e)| {
                rd.extract(n, *s, *e)
                    .map(|b| seqformat::Sequence::new(region_name(n, *s, *e), b))
            })
            .collect::<Result<_>>()?;
        if o.flags.contains("http-stats") {
            if let Some((reqs, bytes)) = rd.http_stats() {
                eprintln!("http: {reqs} requests, {bytes} bytes");
            }
        }
        write_extract_output(o.vals.get("out"), &out_seqs, width)?;
        return Ok(());
    }

    // Peek only a small prefix for format detection so an indexed file isn't
    // slurped into memory just to fetch one region — the reader then opens it
    // seek-based. A few regions (per-fetch / "grab one contig") open seek-based
    // for minimal latency; a large batch or extract-all amortises a single
    // whole-file read, which beats thousands of per-region seeks.
    let prefix = read_prefix(file, 64)?;
    let bulk = regions.is_empty() || regions.len() > 1024;
    let out_seqs: Vec<seqformat::Sequence> = if twobit::is_twobit(&prefix) {
        if bulk {
            run!(twobit::TwoBitReader::from_vec(std::fs::read(file)?)?)
        } else {
            run!(twobit::TwoBitReader::open(file)?)
        }
    } else if twobyte::is_twobyte(&prefix) {
        // 2be carries an on-disk B+ tree, so a seek-open lookup is O(log N).
        if bulk {
            run!(twobyte::TwoByteReader::from_vec(std::fs::read(file)?)?)
        } else {
            run!(twobyte::TwoByteReader::open(file)?)
        }
    } else if fourbit::is_fourbit(&prefix) {
        // 4bit has no index: a lookup must touch every record header either way,
        // so one sequential slurp beats O(N) scattered seeks.
        run!(fourbit::FourBitReader::from_vec(std::fs::read(file)?)?)
    } else if samtools::is_bgzf(&prefix) || samtools::is_fasta(&prefix) {
        run!(samtools::FaidxReader::open(file)?)
    } else {
        run!(twobit::TwoBitReader::from_vec(std::fs::read(file)?)?)
    };

    write_extract_output(o.vals.get("out"), &out_seqs, width)
}

/// Write extracted sequences to `--out` (FASTA file) or stdout.
fn write_extract_output(
    out: Option<&String>,
    seqs: &[seqformat::Sequence],
    width: usize,
) -> Result<()> {
    match out {
        Some(path) => fasta::write_file(path, seqs, width)?,
        None => {
            use std::io::Write;
            let bytes = fasta::write_bytes(seqs, width);
            std::io::stdout().write_all(&bytes)?;
        }
    }
    Ok(())
}

fn cmd_random(args: &[String]) -> Result<()> {
    let o = parse(
        args,
        &[],
        &[
            "seqs", "length", "n-frac", "iub-frac", "seed", "width", "prefix", "n-runs",
            "iub-runs",
        ],
    )?;
    let output = o
        .pos
        .first()
        .ok_or_else(|| Error::Format("random needs an output file".into()))?;
    let opts = GenOpts {
        seqs: o.val_usize("seqs", 1)?,
        length: o.val_usize("length", 1_000_000)?,
        n_frac: o.val_f64("n-frac", 0.01)?,
        iub_frac: o.val_f64("iub-frac", 0.005)?,
        seed: o.val_u64("seed", 1)?,
        name_prefix: o.vals.get("prefix").cloned().unwrap_or_else(|| "seq".into()),
        n_runs: o.val_usize("n-runs", 0)?,
        iub_runs: o.val_usize("iub-runs", 0)?,
    };
    let seqs = generate(&opts);
    fasta::write_file(output, &seqs, o.val_usize("width", 60)?)?;
    let mode = |runs: usize| if runs == 0 { "scattered".to_string() } else { format!("{runs} runs") };
    eprintln!(
        "wrote {} sequence(s) x {} bp to {output} (N {} [{}], IUB {} [{}], seed {})",
        opts.seqs,
        opts.length,
        opts.n_frac,
        mode(opts.n_runs),
        opts.iub_frac,
        mode(opts.iub_runs),
        opts.seed
    );
    Ok(())
}

fn cmd_info(args: &[String]) -> Result<()> {
    let o = parse(args, &[], &[])?;
    let path = o
        .pos
        .first()
        .ok_or_else(|| Error::Format("info needs a file argument".into()))?;
    let data = std::fs::read(path)?;

    if twobyte::is_twobyte(&data) {
        let rd = twobyte::TwoByteReader::from_vec(data)?;
        let stats = rd.sequence_stats()?;
        println!("{path}: 2be (B+ tree TOC + merged tagged-edit stream)");
        println!("  sequences: {}", stats.len());
        println!("  {:<24} {:>12} {:>10}", "name", "length", "edits");
        for (name, len, edits) in stats {
            println!("  {:<24} {:>12} {:>10}", name, len, edits);
        }
    } else if fourbit::is_fourbit(&data) {
        let seqs = fourbit::from_bytes(&data)?;
        println!("{path}: 4-bit sequence file (BWA/BAM nibble encoding)");
        println!("  sequences: {}", seqs.len());
        println!("  {:<24} {:>12}", "name", "length");
        for s in &seqs {
            println!("  {:<24} {:>12}", s.name, s.bases.len());
        }
    } else if samtools::is_bgzf(&data) || samtools::is_fasta(&data) {
        let rd = samtools::FaidxReader::open(path)?;
        let kind = if samtools::is_bgzf(&data) {
            "BGZF-compressed FASTA (samtools .fai/.gzi)"
        } else {
            "FASTA (samtools .fai)"
        };
        println!("{path}: {kind}");
        println!("  sequences: {}", rd.names().len());
        println!("  {:<24} {:>12}", "name", "length");
        for (name, len) in rd.sequence_infos() {
            println!("  {:<24} {:>12}", name, len);
        }
    } else {
        let tb = twobit::from_bytes(&data)?;
        println!(
            "{path}: twoBit {}{}{}",
            if tb.long { "long (v1)" } else { "standard (v0)" },
            if tb.iub { " + IUB extension" } else { "" },
            if tb.indexed { " + name index" } else { "" }
        );
        println!("  sequences: {}", tb.sequences.len());
        println!(
            "  {:<20} {:>12} {:>9} {:>9} {:>9}",
            "name", "length", "N-blocks", "mask", "IUB"
        );
        for st in twobit::stats(&tb) {
            println!(
                "  {:<20} {:>12} {:>9} {:>9} {:>9}",
                st.name, st.len, st.n_blocks, st.mask_blocks, st.iub_blocks
            );
        }
    }
    Ok(())
}

/// Read up to `n` leading bytes of a file (fewer if it is shorter) for cheap
/// format detection, without loading the whole file.
fn read_prefix(path: &str, n: usize) -> Result<Vec<u8>> {
    use std::io::Read;
    let mut f = std::fs::File::open(path)?;
    let mut buf = vec![0u8; n];
    let mut filled = 0;
    while filled < n {
        match f.read(&mut buf[filled..])? {
            0 => break,
            k => filled += k,
        }
    }
    buf.truncate(filled);
    Ok(buf)
}

fn two<'a>(o: &'a Opts, what: &str) -> Result<(&'a str, &'a str)> {
    if o.pos.len() != 2 {
        return fail(format!("{what} expects exactly 2 file arguments"));
    }
    Ok((&o.pos[0], &o.pos[1]))
}
