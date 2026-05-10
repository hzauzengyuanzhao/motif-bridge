# Perl implementation

This directory contains the Perl 5 implementation of `meme2homer` and `homer2meme`.

Use this implementation when you are working on a server where Perl is already available and you do not want to compile Rust binaries or configure a Python environment.

## Commands

```bash
perl perl_scripts/meme2homer.pl -i motifs.meme -j JASPAR2026 > motifs.homer
perl perl_scripts/homer2meme.pl -i motifs.homer > motifs.meme
```

Both scripts support plain text input, `.gz` input, and stdin through `-i -`.

## Notes

The Perl implementation is intended to match the main CLI behavior of the Python and Rust implementations for standard DNA A/C/G/T motif files. See the top-level README for metadata and round-trip limitations.
