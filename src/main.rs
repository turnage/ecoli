extern crate hmm;

mod corpus;

use std::collections::HashMap;
use std::io::Read;
use std::error::Error;

use hmm::base::Emitter;

fn file_string(filename: &str) -> Result<String, std::io::Error> {
    let mut tfile = std::fs::File::open(filename)?;
    let mut buffer = String::new();
    let _ = tfile.read_to_string(&mut buffer)?;
    Ok(buffer)
}

fn sequence_file(filename: &str) -> Result<Vec<corpus::Sequence>, String> {
    file_string(filename)
        .map_err(|e| String::from(e.description()))
        .and_then(|file| {
            file.lines()
                .map(|line| line.parse::<corpus::Sequence>())
                .collect()
        })
}

fn main() {
    let bases = corpus::parse_dna(file_string("resources/ecoli_m54.txt").expect("dna file"))
        .expect("parsed dna bases");

    let train_ranges = sequence_file("resources/train.txt").expect("parsed train sequences");
    let test_ranges = sequence_file("resources/test.txt").expect("parsed test sequences");

    let mut internal_codon_table = HashMap::new();
    let mut stop_codon_table = HashMap::new();
    let train_seqs = train_ranges.iter()
        .map(|r| r.label_dna(&bases, &mut internal_codon_table, &mut stop_codon_table))
        .collect::<Vec<Vec<(corpus::Label, corpus::Base)>>>();
    let tunings = corpus::Sequence::trans_tunings(&internal_codon_table, &stop_codon_table);
    let model = hmm::base::train::discrete(&train_seqs, None);

    let ground_truth = test_ranges.iter()
        .map(|r| r.label_dna(&bases, &mut internal_codon_table, &mut stop_codon_table))
        .collect::<Vec<Vec<(corpus::Label, corpus::Base)>>>();
    let test_seqs = test_ranges.iter().map(|r| r.dna(&bases)).collect::<Vec<Vec<corpus::Base>>>();
    let results = test_seqs.iter()
    //    .take(1)
        .map(|s| hmm::base::Solve::most_probable_sequence(s, &model))
        .collect::<Result<Vec<Vec<corpus::Label>>, String>>()
        .expect("prediction results");
    // for (i, result) in results.iter().enumerate() {
    // for (j, label) in result.iter().zip(ground_truth[i].iter()).enumerate() {
    // println!("{}: {:?}", j, label);
    // }
    // }

    let coords = results.iter()
        .enumerate()
        .map(|(i, labels)| corpus::Sequence::from_labels(test_ranges[i].start(), labels))
        .collect::<Result<Vec<corpus::Sequence>, String>>()
        .expect("sequence parsed prediction results");
    for coord in coords.iter() {
        println!("{}", coord);
    }
}
