#[allow(unused_imports)]
use super::prelude::*;

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Nucleotide {
    A,
    C,
    G,
    T,
    Unknown,
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
            'n' => Unknown,
            _ => panic!("Invalid character : {:?}", ch),
        }
    }
}

impl Into<char> for Nucleotide {
    fn into(self) -> char {
        use Nucleotide::*;
        match self {
            A => 'a',
            C => 'c',
            G => 'g',
            T => 't',
            Unknown => 'n',
        }
    }
}
