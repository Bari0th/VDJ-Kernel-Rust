#[allow(unused_imports)]
use super::prelude::*;
use std::cmp::PartialEq;

#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq)]
pub enum Nucleotide {
    A = 'a' as isize,
    C = 'c' as isize,
    G = 'g' as isize,
    T = 't' as isize,
    N = 'n' as isize,
}

use std::convert::From;

impl From<char> for Nucleotide {
    fn from(ch: char) -> Nucleotide {
        use Nucleotide::*;
        match ch {
            'A' | 'a' => A,
            'C' | 'c' => C,
            'G' | 'g' => G,
            'T' | 't' => T,
            'n' | 'N' => N,
            _ => panic!("Invalid character : {:?}", ch),
        }
    }
}
