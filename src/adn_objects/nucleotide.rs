#[allow(unused_imports)]
use super::prelude::*;
use std::cmp::PartialEq;

#[derive(Debug, Clone, Copy)]
pub enum Nucleotide {
    A = 'a' as isize,
    C = 'c' as isize,
    G = 'g' as isize,
    T = 't' as isize,
    AC = 0,
    AG,
    AT,
    CG,
    CT,
    GT,
    ACG,
    ACT,
    AGT,
    CGT,
    ACGT = 'n' as isize,
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
            'n' | 'N' => ACGT,
            _ => panic!("Invalid character : {:?}", ch),
        }
    }
}

impl PartialEq for Nucleotide {
    fn eq(&self, other: &Nucleotide) -> bool {
        use Nucleotide::*;
        match other {
            A | C | G | T => match self {
                A | C | G | T => *self as isize == *other as isize,
                AC => *other as isize == 'c' as isize || *other as isize == 'a' as isize,
                AG => *other as isize == 'g' as isize || *other as isize == 'a' as isize,
                AT => *other as isize == 't' as isize || *other as isize == 'a' as isize,
                CG => *other as isize == 'g' as isize || *other as isize == 'a' as isize,
                CT => *other as isize == 't' as isize || *other as isize == 'a' as isize,
                GT => *other as isize == 't' as isize || *other as isize == 'a' as isize,
                ACG => *other as isize != 't' as isize,
                ACT => *other as isize != 'g' as isize,
                AGT => *other as isize != 'c' as isize,
                CGT => *other as isize != 'a' as isize,
                ACGT => true,
            },

            ACGT => true,

            _ => false,
        }
    }
}
