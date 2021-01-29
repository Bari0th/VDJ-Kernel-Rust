use super::prelude::*;
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

    pub fn increment(&mut self, prediction: &str, result: &str) {
        if let Some(&idx) = self.1.get(prediction) {
            let i = idx;
            if let Some(&idx) = self.1.get(result) {
                let j = idx;
                self.0[i][j] += 1;
            }
        }
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
