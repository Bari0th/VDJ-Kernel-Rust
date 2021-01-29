#![allow(dead_code)]

use test_some_algorithm::adn_objects::prelude::*;

fn main() {
    let genes_j = Genes::load_from("data/database/IGHJgenes.txt");
    let genes_d = Genes::load_from("data/database/IGHDgenes.txt");
    let genes_v = Genes::load_from("data/database/IGHVgenes.txt");

    let sequences = Sequences::load_from("data/datasets/simulations/sim1.fasta");

    let results = SequenceResult::load_from("data/datasets/simulations/simTrueGenes.txt");

    let predictions = SequenceResult::calcul_and_save_to(
        "test_alignment/genes_alignment.txt",
        &sequences,
        &genes_v,
        &genes_d,
        &genes_j,
        &results,
    );
    // let predictions = SequenceResult::load_from("../data/datasets/simulations/simTrueGenes.txt");

    predictions.save_to("first_predictions");

    let confusion_matrixes = SequenceResult::get_confusion_matrixes(
        &predictions,
        &results,
        &genes_v,
        &genes_d,
        &genes_j,
    );

    println!(
        "{}\n\n{}\n\n{}",
        confusion_matrixes.0, confusion_matrixes.1, confusion_matrixes.2
    );
}
