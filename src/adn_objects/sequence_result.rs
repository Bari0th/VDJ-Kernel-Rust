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
                .map(|(key, value)| {
                    format!(
                        "{:8}\t{:20}\t{:20}\t{:20}\n",
                        key, value.0, value.1, value.2
                    )
                })
                .collect::<String>()
                .as_bytes(),
        )
        .expect("Cannot write to file");
    }

    pub fn compare_with(&self, results: &SequenceResult, path: &str) {
        use std::fs::File;
        use std::io::BufWriter;
        use std::io::Write;

        let mut writer = BufWriter::new(File::create(path).expect("Cannot create file"));

        let mut keys = self.0.keys().collect::<Vec<_>>();
        keys.sort_by_key(|key| key.strip_prefix("S").unwrap().parse::<usize>().unwrap());

        for key in keys {
            if let (Some(res), Some(pred)) = (results.0.get(key), self.0.get(key)) {
                write!(
                    writer,
                    "{:8}Pred : {:20}{:20}{:20}\n        Res  : {:20}{:20}{:20}\n\n",
                    key, pred.0, pred.1, pred.2, res.0, res.1, res.2
                )
                .expect("Cannot write to file");
            }
        }
    }

    pub fn confusion_with(&self, results: &SequenceResult, path: &str) {
        use std::fs::File;
        use std::io::Write;

        let mut file = File::create(path).expect("Cannot create file");

        let (v, d, j, n) =
            self.0
                .iter()
                .fold((0, 0, 0, 0), |(mut v, mut d, mut j, mut n), (key, pred)| {
                    if let Some(res) = results.0.get(key) {
                        if res.0 == pred.0 {
                            v += 1;
                        }
                        if res.1 == pred.1 {
                            d += 1;
                        }
                        if res.2 == pred.2 {
                            j += 1;
                        }
                        n += 1;
                    }
                    (v, d, j, n)
                });

        let v_per = v as f64 / n as f64 * 100.;
        let d_per = d as f64 / n as f64 * 100.;
        let j_per = j as f64 / n as f64 * 100.;
        let tot = v + d + j;
        let tot_per = tot as f64 / (3 * n) as f64 * 100.;

        writeln!(file, "V  P : {:4} | {:3.2}%", v, v_per).expect("Cannot write to file");
        writeln!(file, "   F : {:4} | {:3.2}%\n", n - v, 100. - v_per)
            .expect("Cannot write to file");
        writeln!(file, "D  P : {:4} | {:3.2}%", d, d_per).expect("Cannot write to file");
        writeln!(file, "   F : {:4} | {:3.2}%\n", n - d, 100. - d_per)
            .expect("Cannot write to file");
        writeln!(file, "J  P : {:4} | {:3.2}%", j, j_per).expect("Cannot write to file");
        writeln!(file, "   F : {:4} | {:3.2}%\n", n - j, 100. - j_per)
            .expect("Cannot write to file");
        writeln!(file, "ALL P : {:4} | {:3.2}%", tot, tot_per).expect("Cannot write to file");
        writeln!(
            file,
            "    F : {:4} | {:3.2}%\n",
            3 * n - tot,
            100. - tot_per
        )
        .expect("Cannot write to file");
    }

    pub fn confusion_with_without_allele(&self, results: &SequenceResult, path: &str) {
        use std::fs::File;
        use std::io::Write;

        let mut file = File::create(path).expect("Cannot create file");

        let (v, d, j, n) =
            self.0
                .iter()
                .fold((0, 0, 0, 0), |(mut v, mut d, mut j, mut n), (key, pred)| {
                    if let Some(res) = results.0.get(key) {
                        if res.0.split('*').next() == pred.0.split('*').next() {
                            v += 1;
                        }
                        if res.1.split('*').next() == pred.1.split('*').next() {
                            d += 1;
                        }
                        if res.2.split('*').next() == pred.2.split('*').next() {
                            j += 1;
                        }
                        n += 1;
                    }
                    (v, d, j, n)
                });
        let v_per = v as f64 / n as f64 * 100.;
        let d_per = d as f64 / n as f64 * 100.;
        let j_per = j as f64 / n as f64 * 100.;
        let tot = v + d + j;
        let tot_per = tot as f64 / (3 * n) as f64 * 100.;
        writeln!(file, "V  P : {:4} | {:3.2}%", v, v_per).expect("Cannot write to file");
        writeln!(file, "   F : {:4} | {:3.2}%\n", n - v, 100. - v_per)
            .expect("Cannot write to file");
        writeln!(file, "D  P : {:4} | {:3.2}%", d, d_per).expect("Cannot write to file");
        writeln!(file, "   F : {:4} | {:3.2}%\n", n - d, 100. - d_per)
            .expect("Cannot write to file");
        writeln!(file, "J  P : {:4} | {:3.2}%", j, j_per).expect("Cannot write to file");
        writeln!(file, "   F : {:4} | {:3.2}%\n", n - j, 100. - j_per)
            .expect("Cannot write to file");
        writeln!(file, "ALL P : {:4} | {:3.2}%", tot, tot_per).expect("Cannot write to file");
        writeln!(
            file,
            "    F : {:4} | {:3.2}%\n",
            3 * n - tot,
            100. - tot_per
        )
        .expect("Cannot write to file");
    }

    #[allow(clippy::type_complexity)]
    pub fn get_best_genes<'a>(
        sequence: &(String, Sequence),
        genes_v: &'a Genes,
        genes_d: &'a Genes,
        genes_j: &'a Genes,
    ) -> ((&'a str, usize), (&'a str, usize), (&'a str, usize)) {
        let best_v = genes_v.get_best::<10>(&sequence.1);
        let best_d = genes_d.get_best_all_sequence(&sequence.1[best_v.1 - 10..]);
        let best_j = genes_j.get_best_all_sequence(&sequence.1[best_v.1 + best_d.1 - 20..]);
        (best_v, best_d, best_j)
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
                    //println!("{}", sequence.0);
                    let (best_v, best_d, best_j) =
                        SequenceResult::get_best_genes(sequence, genes_v, genes_d, genes_j);
                    (
                        sequence.0.clone(),
                        (
                            best_v.0.to_string(),
                            best_d.0.to_string(),
                            best_j.0.to_string(),
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
                    let (best_v, best_d, best_j) =
                        SequenceResult::get_best_genes(sequence, genes_v, genes_d, genes_j);
                    let tmp = (
                        sequence.0.clone(),
                        (
                            best_v.0.to_string(),
                            best_d.0.to_string(),
                            best_j.0.to_string(),
                        ),
                    );
                    let resultat = resultats.0.get(&tmp.0).unwrap();
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
                    let (best_v, best_d, best_j) =
                        SequenceResult::get_best_genes(sequence, genes_v, genes_d, genes_j);
                    let tmp = (
                        sequence.0.clone(),
                        (
                            best_v.0.to_string(),
                            best_d.0.to_string(),
                            best_j.0.to_string(),
                        ),
                    );
                    println!("{}", sequence.0);
                    let result = results.0.get(&tmp.0).unwrap();

                    let seq = sequence.1.to_string();
                    let seq_bytes = seq.as_bytes();

                    let v_pred = Sequence::get_best_alignment(
                        genes_v.get_gene(&tmp.1 .0).unwrap(),
                        &sequence.1,
                    )
                    .1;

                    let d_pred = ('.'..='.').cycle().take(best_v.1 - 10).collect::<String>()
                        + &Sequence::get_best_alignment(
                            genes_d.get_gene(&tmp.1 .1).unwrap(),
                            &sequence.1[best_v.1 - 10..],
                        )
                        .1;

                    let j_pred = ('.'..='.')
                        .cycle()
                        .take(best_v.1 + best_d.1 - 20)
                        .collect::<String>()
                        + &Sequence::get_best_alignment(
                            genes_j.get_gene(&tmp.1 .2).unwrap(),
                            &sequence.1[best_v.1 + best_d.1 - 20..],
                        )
                        .1;

                    let v_res = if let Some(gene) = genes_v.get_gene(&result.0) {
                        Sequence::get_best_alignment(gene, &sequence.1).1
                    } else {
                        ('-'..='-')
                            .cycle()
                            .take(sequence.1.len())
                            .collect::<String>()
                    };

                    let d_res = if let Some(gene) = genes_d.get_gene(&result.1) {
                        Sequence::get_best_alignment(gene, &sequence.1).1
                    } else {
                        ('-'..='-')
                            .cycle()
                            .take(sequence.1.len())
                            .collect::<String>()
                    };

                    let j_res = if let Some(gene) = genes_j.get_gene(&result.2) {
                        Sequence::get_best_alignment(gene, &sequence.1).1
                    } else {
                        ('-'..='-')
                            .cycle()
                            .take(sequence.1.len())
                            .collect::<String>()
                    };

                    let preds: [&str; 3] = [&v_pred, &d_pred, &j_pred];
                    let res: [&str; 3] = [&v_res, &d_res, &j_res];

                    let check = |i: usize, seqs: [&str; 3]| {
                        Ok(true)
                            == seqs
                                .iter()
                                .map(|seq| seq.as_bytes())
                                .try_fold(false, |b, s| match s[i] {
                                    b'.' => Ok(b),
                                    b'-' => Ok(true),
                                    c if c == seq_bytes[i] => {
                                        if !b {
                                            Ok(true)
                                        } else {
                                            Err(())
                                        }
                                    }
                                    _ => Err(()),
                                })
                    };

                    let indices: Vec<_> = (0..seq_bytes.len())
                        .filter(|&i| !check(i, preds) || !check(i, res))
                        .collect();

                    let seq_to_string = |idx: &[usize], seq: &str| -> String {
                        let mut s = if idx[0] > 6 {
                            format!("[{}]-{}", idx[0], &seq[idx[0]..=idx[0]])
                        } else {
                            (&seq[0..=idx[0]]).to_string()
                        };

                        let mut last_i = idx[0];

                        for &i in idx.iter().skip(1) {
                            if i > last_i + 7 {
                                s += &format!("-[{}]-{}", i - last_i - 1, &seq[i..=i]);
                            } else {
                                s += &seq[last_i + 1..=i];
                            }
                            last_i = i;
                        }

                        if seq.len() - 1 > last_i + 6 {
                            s += &format!("-[{}]\n", seq.len() - last_i - 2);
                        } else {
                            s = s + &seq[last_i + 1..] + "\n";
                        };
                        s
                    };

                    if indices.is_empty() {
                        file.write_all(
                            format!("Sequence : {}\n--Matches Perfectly--\n\n", tmp.0).as_bytes(),
                        )
                        .expect("Cannot write to file");
                    } else {
                        file.write_all(
                            format!(
                                "Sequence : {}\n{}\nPredictions :\n{}\n\nResults :\n{}\n\n\n",
                                tmp.0,
                                seq_to_string(&indices, &seq),
                                preds
                                    .iter()
                                    .map(|seq| seq_to_string(&indices, seq))
                                    .collect::<String>(),
                                res.iter()
                                    .map(|seq| seq_to_string(&indices, seq))
                                    .collect::<String>(),
                            )
                            .as_bytes(),
                        )
                        .expect("Cannot write to file");
                    }
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

    pub fn get_confusion_matrixes_without_allele(
        predictions: &SequenceResult,
        results: &SequenceResult,
        genes_v: &Genes,
        genes_d: &Genes,
        genes_j: &Genes,
    ) -> (ConfusionMatrix, ConfusionMatrix, ConfusionMatrix) {
        let mut confusion_v = ConfusionMatrix::new_without_allele(genes_v);
        let mut confusion_d = ConfusionMatrix::new_without_allele(genes_d);
        let mut confusion_j = ConfusionMatrix::new_without_allele(genes_j);
        predictions.0.iter().for_each(|(name, prediction)| {
            let result = results
                .0
                .get(name)
                .expect("Cannot get result for this sequence");
            confusion_v.increment_without_allele(&prediction.0, &result.0);
            confusion_d.increment_without_allele(&prediction.1, &result.1);
            confusion_j.increment_without_allele(&prediction.2, &result.2);
        });
        (confusion_v, confusion_d, confusion_j)
    }

    pub fn calcul_with_k_mers<
        'a,
        const KV: usize,
        const KD: usize,
        const KJ: usize,
        const RANGEV: isize,
        const RANGED: isize,
        const RANGEJ: isize,
    >(
        sequences: &Sequences,
        genes_v: &'a Genes,
        genes_d: &'a Genes,
        genes_j: &'a Genes,
    ) -> SequenceResult {
        let v_tree = genes_v.construct_k_mers_trees::<KV>();
        let d_tree = genes_d.construct_k_mers_trees::<KD>();
        let j_tree = genes_j.construct_k_mers_trees::<KJ>();
        SequenceResult(
            sequences
                .iter()
                .map(|sequence| {
                    let best_v = sequence
                        .1
                        .slice(0..200.min(sequence.1.len()))
                        .get_best_matches::<KV, RANGEV>(&v_tree)[0]
                        .0;
                    let size_v = genes_v.get_gene(best_v).unwrap().len();

                    let best_j = sequence
                        .1
                        .slice(size_v + 10..)
                        .get_best_matches::<KJ, RANGEJ>(&j_tree)[0]
                        .0;
                    let size_j = genes_j.get_gene(best_j).unwrap().len();

                    let best_d = sequence
                        .1
                        .slice(size_v..size_v.max(sequence.1.len() - size_j))
                        .get_best_matches::<KD, RANGED>(&d_tree)[0]
                        .0;
                    (
                        sequence.0.clone(),
                        (best_v.to_string(), best_d.to_string(), best_j.to_string()),
                    )
                })
                .collect(),
        )
    }
}
