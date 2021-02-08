use super::prelude::*;
use std::collections::BTreeMap;
use std::collections::HashMap;
use std::ops::Deref;

#[derive(Debug)]
pub struct Genes(Vec<(String, Sequence)>);

impl Genes {
    pub fn get_best<'a, const N: usize>(&'a self, sequence: &[Nucleotide]) -> (&'a str, usize) {
        self.0
            .iter()
            .fold((("No genes", 0), i32::MIN), |res, gene| {
                let score = Sequence::compare::<N>(&gene.1, sequence);
                if score > res.1 {
                    ((&gene.0, gene.1.len()), score)
                } else {
                    res
                }
            })
            .0
    }

    pub fn get_best_all_sequence<'a>(&'a self, sequence: &[Nucleotide]) -> (&'a str, usize) {
        self.0
            .iter()
            .fold((("No genes", 0), i32::MIN), |res, gene| {
                let score = Sequence::compare_all_sequence(&gene.1, sequence);
                if score > res.1 {
                    ((&gene.0, gene.1.len()), score)
                } else {
                    res
                }
            })
            .0
    }

    pub fn get_gene<'a>(&'a self, name: &str) -> Option<&'a Sequence> {
        self.iter()
            .filter_map(|(s, seq)| if s == name { Some(seq) } else { None })
            .next()
    }

    pub fn load_from(path: &str) -> Self {
        use std::fs::File;
        use std::io::Read;
        let mut file = File::open(path).expect("Cannot open file");

        let mut buffer = String::new();

        file.read_to_string(&mut buffer).expect("Cannot read file");

        Genes(
            buffer
                .split('>')
                .skip(1)
                .map(|s| {
                    let mut format_gene = s.split('|');
                    (
                        format_gene.nth(1).unwrap().to_string(),
                        Sequence::new(format_gene.last().unwrap()),
                    )
                })
                .collect(),
        )
    }

    pub fn construct_k_mers_trees<const K: usize>(
        &self,
    ) -> HashMap<&[Nucleotide], BTreeMap<usize, Vec<&str>>> {
        let mut k_mers_trees: HashMap<&[Nucleotide], BTreeMap<usize, Vec<&str>>> = HashMap::new();

        for (name, sequence) in self.0.iter() {
            let k_mers = sequence.decompose_k_mers::<K>();
            for (k_mer, &count) in k_mers.iter() {
                if let Some(k_mers_tree) = k_mers_trees.get_mut(k_mer) {
                    if let Some(names) = k_mers_tree.get_mut(&count) {
                        names.push(name);
                    } else {
                        k_mers_tree.insert(count, vec![name.as_ref()]);
                    }
                } else {
                    let mut k_mers_tree = BTreeMap::new();
                    k_mers_tree.insert(count, vec![name.as_ref()]);
                    k_mers_trees.insert(k_mer, k_mers_tree);
                }
            }
        }

        k_mers_trees
    }
}

impl Deref for Genes {
    type Target = Vec<(String, Sequence)>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
