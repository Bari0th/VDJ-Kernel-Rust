#![allow(dead_code)]

use std::time::Instant;
use test_some_algorithm::adn_objects::prelude::*;

fn main() {
    let genes_j = Genes::load_from("data/database/IGHJgenes.txt");
    let genes_d = Genes::load_from("data/database/IGHDgenes.txt");
    let genes_v = Genes::load_from("data/database/IGHVgenes.txt");

    let sequences = Sequences::load_from("data/datasets/simulations/sim1.fasta");

    let results = SequenceResult::load_from("data/datasets/simulations/simTrueGenes.txt");

    let start = Instant::now();
    let predictions = SequenceResult::calcul_from(&sequences, &genes_v, &genes_d, &genes_j);
    // let predictions = SequenceResult::calcul_and_save_to(
    //     "test_alignment/genes_alignment.txt",
    //     &sequences,
    //     &genes_v,
    //     &genes_d,
    //     &genes_j,
    //     &results,
    // );
    let time = start.elapsed().as_secs_f64();
    //predictions.save_to("my_prediction.txt");

    println!(
        "Time for 1 000 000 sequences : {}",
        time / sequences.len() as f64 * 1_000_000.
    );
    //let predictions = SequenceResult::load_from("my_predictions_tmp.txt");

    predictions.save_to("my_predictions_tmp.txt");

    predictions.compare_with(&results, "result_comparison.txt");

    predictions.confusion_with(&results, "confusion_matrix.txt");

    predictions.confusion_with_without_allele(&results, "confusion_matrix_without_allele.txt");

    let confusion_matrixes = SequenceResult::get_confusion_matrixes_without_allele(
        &predictions,
        &results,
        &genes_v,
        &genes_d,
        &genes_j,
    );
    let conf_v = DiffList::from_conf(&confusion_matrixes.0);
    let conf_d = DiffList::from_conf(&confusion_matrixes.1);
    let conf_j = DiffList::from_conf(&confusion_matrixes.2);

    // println!("{}\n\n{}\n\n{}", conf_v, conf_d, conf_j);

    // println!("{:?}\n\n{:?}\n\n{:?}", conf_v, conf_d, conf_j);

    use std::fs::File;
    use std::io::Write;

    File::create("accuracy.txt")
        .expect("Cannot create file")
        .write_all(format!("{}\n\n{}\n\n{}", conf_v, conf_d, conf_j).as_bytes())
        .expect("Cannot write to file");

    File::create("accuracy_most_choose.txt")
        .expect("Cannot create file")
        .write_all(format!("{:?}\n\n{:?}\n\n{:?}", conf_v, conf_d, conf_j).as_bytes())
        .expect("Cannot write to file");
}
