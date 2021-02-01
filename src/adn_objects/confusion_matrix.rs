use super::prelude::*;
use core::ops::Deref;
use std::collections::HashMap;
use std::fmt;

pub struct ConfusionMatrix(Vec<Vec<usize>>, HashMap<String, usize>);

impl ConfusionMatrix {
    pub fn new(genes: &Genes) -> ConfusionMatrix {
        let (size, hash_index) =
            genes
                .iter()
                .fold((0, HashMap::new()), |(idx, mut hash_map), gene| {
                    if hash_map.contains_key(&gene.0) {
                        (idx, hash_map)
                    } else {
                        hash_map.insert(gene.0.clone(), idx);
                        (idx + 1, hash_map)
                    }
                });
        let confusion = vec![vec![0; size]; size];
        ConfusionMatrix(confusion, hash_index)
    }
    pub fn new_without_allele(genes: &Genes) -> ConfusionMatrix {
        let (size, hash_index) =
            genes
                .iter()
                .fold((0, HashMap::new()), |(idx, mut hash_map), gene| {
                    if hash_map.contains_key(gene.0.split('*').next().unwrap()) {
                        (idx, hash_map)
                    } else {
                        hash_map.insert(gene.0.split('*').next().unwrap().to_string(), idx);
                        (idx + 1, hash_map)
                    }
                });
        let confusion = vec![vec![0; size]; size];
        ConfusionMatrix(confusion, hash_index)
    }

    pub fn get_ids(&self) -> Vec<String> {
        let mut res = self
            .1
            .iter()
            .map(|(k, &v)| (v, k.clone()))
            .collect::<Vec<_>>();
        res.sort();
        res.into_iter().map(|(_, k)| k).collect()
    }

    pub fn increment(&mut self, prediction: &str, result: &str) {
        if let Some(&idx) = self.1.get(prediction) {
            let i = idx;
            if let Some(&idx) = self.1.get(result) {
                let j = idx;
                self.0[i][j] += 1;
            }
        }
    }

    pub fn increment_without_allele(&mut self, prediction: &str, result: &str) {
        if let Some(&idx) = self.1.get(prediction.split('*').next().unwrap()) {
            let i = idx;
            if let Some(&idx) = self.1.get(result.split('*').next().unwrap()) {
                let j = idx;
                self.0[i][j] += 1;
            }
        }
    }
}

impl Deref for ConfusionMatrix {
    type Target = [Vec<usize>];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for ConfusionMatrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            self.0
                .iter()
                .map(|line| format!("{:?}\n", line))
                .collect::<String>()
        )
    }
}
