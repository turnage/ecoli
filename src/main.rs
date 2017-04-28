#[macro_use]
extern crate clap;
extern crate hmm;

mod corpus;

use std::collections::HashMap;
use std::io::Read;
use std::error::Error;

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

fn deinstance_label(la: corpus::Label) -> corpus::Label {
    match la {
        corpus::Label::StopCodon1(_) => corpus::Label::StopCodon1(0),
        corpus::Label::StopCodon2(_) => corpus::Label::StopCodon2(0),
        corpus::Label::StopCodon3(_) => corpus::Label::StopCodon3(0),
        corpus::Label::InternalCodon1(_) => corpus::Label::InternalCodon1(0),
        corpus::Label::InternalCodon2(_) => corpus::Label::InternalCodon2(0),
        corpus::Label::InternalCodon3(_) => corpus::Label::InternalCodon3(0),
        _ => la,
    }
}

fn main() {
    let yaml = load_yaml!("cli.yaml");
    let matches = clap::App::from_yaml(yaml).get_matches();

    let dna_filename = matches.value_of("dna").unwrap();
    let truth_filename = matches.value_of("truth").unwrap();
    let seq_filename = matches.value_of("sequences").unwrap();

    let bases = corpus::parse_dna(file_string(&dna_filename).expect("dna file"))
        .expect("parsed dna bases");

    let train_ranges = sequence_file(&truth_filename).expect("parsed train sequences");
    let test_ranges = sequence_file(&seq_filename).expect("parsed test sequences");

    let mut internal_codon_table = HashMap::new();
    let mut stop_codon_table = HashMap::new();
    let train_seqs = train_ranges.iter()
        .map(|r| r.label_dna(&bases, &mut internal_codon_table, &mut stop_codon_table))
        .collect::<Vec<Vec<(corpus::Label, corpus::Base)>>>();
    let model = hmm::base::train::discrete(&train_seqs, None);

    let mut trans_count = HashMap::new();
    let mut out_count = HashMap::new();
    for (&s1, dist) in model.trans.iter() {
        let from = deinstance_label(s1);
        for (&s2, p) in dist.iter() {
            let to = deinstance_label(s2);
            *out_count.entry(from).or_insert(0f64) += *p;
            *trans_count.entry((from, to)).or_insert(0f64) += *p;
        }
    }

    println!("Transition probabilities from training set: ");
    for (&(from, to), &p) in trans_count.iter() {
        println!("{:?} -> {:?}: {}", from, to, p / out_count[&from])
    }

    let test_seqs = test_ranges.iter().map(|r| r.dna(&bases)).collect::<Vec<Vec<corpus::Base>>>();
    let results = test_seqs.iter()
        .map(|s| hmm::base::Solve::most_probable_sequence(s, &model))
        .collect::<Result<Vec<Vec<corpus::Label>>, String>>()
        .expect("prediction results");

    let coords = results.iter()
        .enumerate()
        .map(|(i, labels)| corpus::Sequence::from_labels(test_ranges[i].start(), labels))
        .collect::<Result<Vec<corpus::Sequence>, String>>()
        .expect("sequence parsed prediction results");
    for coord in coords.iter() {
        println!("{}", coord);
    }
}
