use super::prelude::*;
use std::collections::HashMap;

#[derive(Debug)]
pub struct SequenceResult(HashMap<String, (String, String, String)>);

impl SequenceResult {
    pub fn load_from(path: &str) -> Self {
        use std::fs::File;
        use std::io::Read;
        let mut file = File::open(path).expect("Cannot open file");

        let mut buffer = String::new();

        file.read_to_string(&mut buffer).expect("Cannot read file");

        SequenceResult(
            buffer
                .split('\n')
                .rev()
                .skip(1)
                .map(|s| {
                    let mut format_result = s.trim().split('\t');
                    (
                        format_result
                            .next()
                            .expect("Problem reading formated date !")
                            .to_string(),
                        (
                            format_result
                                .next()
                                .expect("Problem reading formated date !")
                                .to_string(),
                            format_result
                                .next()
                                .expect("Problem reading formated date !")
                                .to_string(),
                            format_result
                                .next()
                                .expect("Problem reading formated date !")
                                .to_string(),
                        ),
                    )
                })
                .collect(),
        )
    }

    pub fn save_to(&self, path: &str) {
        use std::fs::File;
        use std::io::Write;
        let mut file = File::create(path).expect("Cannot create file");
        file.write_all(
            self.0
                .iter()
                .map(|(key, value)| format!("{}\t{}\t{}\t{}\n", key, value.0, value.1, value.2))
                .collect::<String>()
                .as_bytes(),
        )
        .expect("Cannot write to file");
    }

    pub fn calcul_from(
        sequences: &Sequences,
        genes_v: &Genes,
        genes_d: &Genes,
        genes_j: &Genes,
    ) -> Self {
        SequenceResult(
            sequences
                .iter()
                .map(|sequence| {
                    (
                        sequence.0.clone(),
                        (
                            genes_v.get_best(&sequence.1).to_string(),
                            genes_d.get_best(&sequence.1).to_string(),
                            genes_j.get_best(&sequence.1).to_string(),
                        ),
                    )
                })
                .collect(),
        )
    }

    pub fn calcul_from_and_compare(
        sequences: &Sequences,
        genes_v: &Genes,
        genes_d: &Genes,
        genes_j: &Genes,
        resultats: &SequenceResult,
    ) -> Self {
        SequenceResult(
            sequences
                .iter()
                .map(|sequence| {
                    let tmp = (
                        sequence.0.clone(),
                        (
                            genes_v.get_best(&sequence.1).to_string(),
                            genes_d.get_best(&sequence.1).to_string(),
                            genes_j.get_best(&sequence.1).to_string(),
                        ),
                    );
                    let resultat = resultats.0.get(&tmp.0.to_string()).unwrap();
                    println!(
                        "{} :\t{} {} {}\n\t{} {} {}\n\n",
                        tmp.0, resultat.0, resultat.1, resultat.2, tmp.1 .0, tmp.1 .1, tmp.1 .2
                    );
                    tmp
                })
                .collect(),
        )
    }

    pub fn calcul_and_save_to(
        path: &str,
        sequences: &Sequences,
        genes_v: &Genes,
        genes_d: &Genes,
        genes_j: &Genes,
        results: &SequenceResult,
    ) -> Self {
        use std::fs::File;
        use std::io::Write;
        let mut file = File::create(path).expect("Cannot create file");
        SequenceResult(
            sequences
                .iter()
                .map(|sequence| {
                    let tmp = (
                        sequence.0.clone(),
                        (
                            genes_v.get_best(&sequence.1).to_string(),
                            genes_d.get_best(&sequence.1).to_string(),
                            genes_j.get_best(&sequence.1).to_string(),
                        ),
                    );
                    println!("{}", sequence.0);
                    let result = results.0.get(&tmp.0.to_string()).unwrap();
                    file.write_all(
                        (format!(
                            "{} :\n{}\n\n{}\n{}\n{}\n\n",
                            tmp.0,
                            sequence.1.to_string(),
                            Sequence::get_best_alignment(
                                genes_v.get_gene(&tmp.1 .0).unwrap(),
                                &sequence.1
                            )
                            .1,
                            Sequence::get_best_alignment(
                                genes_d.get_gene(&tmp.1 .1).unwrap(),
                                &sequence.1
                            )
                            .1,
                            Sequence::get_best_alignment(
                                genes_j.get_gene(&tmp.1 .2).unwrap(),
                                &sequence.1
                            )
                            .1,
                        ) + &format!(
                            "{}\n{}\n{}\n\n",
                            if let Some(gene) = genes_v.get_gene(&result.0) {
                                Sequence::get_best_alignment(gene, &sequence.1).1
                            } else {
                                ('.'..='.')
                                    .cycle()
                                    .take(sequence.1.len())
                                    .collect::<String>()
                            },
                            if let Some(gene) = genes_d.get_gene(&result.1) {
                                Sequence::get_best_alignment(gene, &sequence.1).1
                            } else {
                                ('.'..='.')
                                    .cycle()
                                    .take(sequence.1.len())
                                    .collect::<String>()
                            },
                            if let Some(gene) = genes_j.get_gene(&result.2) {
                                Sequence::get_best_alignment(gene, &sequence.1).1
                            } else {
                                ('.'..='.')
                                    .cycle()
                                    .take(sequence.1.len())
                                    .collect::<String>()
                            },
                        ))
                            .as_bytes(),
                    )
                    .expect("Cannot write to file");
                    tmp
                })
                .collect(),
        )
    }

    pub fn get_confusion_matrixes(
        predictions: &SequenceResult,
        results: &SequenceResult,
        genes_v: &Genes,
        genes_d: &Genes,
        genes_j: &Genes,
    ) -> (ConfusionMatrix, ConfusionMatrix, ConfusionMatrix) {
        let mut confusion_v = ConfusionMatrix::new(genes_v);
        let mut confusion_d = ConfusionMatrix::new(genes_d);
        let mut confusion_j = ConfusionMatrix::new(genes_j);
        predictions.0.iter().for_each(|(name, prediction)| {
            let result = results
                .0
                .get(name)
                .expect("Cannot get result for this sequence");
            confusion_v.increment(&prediction.0, &result.0);
            confusion_d.increment(&prediction.1, &result.1);
            confusion_j.increment(&prediction.2, &result.2);
        });
        (confusion_v, confusion_d, confusion_j)
    }
}
