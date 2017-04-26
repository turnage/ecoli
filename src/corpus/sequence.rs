use std::collections::HashMap;
use std::cmp::Ordering;
use std::fmt;
use std::str::FromStr;

use corpus::base::Base;
use corpus::labeler::Label;
use corpus::labeler;

#[derive(Debug, PartialEq)]
pub struct Sequence {
    start: usize,
    end: usize,
    genes: Vec<(usize, usize)>,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum LabelClass {
    NonCoding,
    StartCodon,
    InternalCodon,
    StopCodon,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
struct LabeledRange {
    class: LabelClass,
    range: (usize, usize),
}

impl PartialOrd for LabeledRange {
    fn partial_cmp(&self, other: &LabeledRange) -> Option<Ordering> {
        self.range.0.partial_cmp(&other.range.0)
    }

    fn lt(&self, other: &LabeledRange) -> bool {
        self.range.0 < other.range.0
    }

    fn le(&self, other: &LabeledRange) -> bool {
        self.range.0 <= other.range.0
    }

    fn gt(&self, other: &LabeledRange) -> bool {
        self.range.0 > other.range.0
    }

    fn ge(&self, other: &LabeledRange) -> bool {
        self.range.0 >= other.range.0
    }
}

impl Ord for LabeledRange {
    fn cmp(&self, other: &LabeledRange) -> Ordering {
        self.range.0.cmp(&other.range.0)
    }
}

impl Sequence {
    pub fn start(&self) -> usize {
        self.start
    }

    pub fn dna(&self, dna: &[Base]) -> Vec<Base> {
        dna[self.start..self.end].to_vec()
    }

    pub fn label_dna(&self,
                     dna: &[Base],
                     internal_codon_table: &mut HashMap<(Base, Base, Base), usize>,
                     stop_codon_table: &mut HashMap<(Base, Base, Base), usize>)
                     -> Vec<(Label, Base)> {
        let mut emissions = Vec::new();
        for range in self.labeled_ranges().iter() {
            let (start, end) = range.range;
            let bases = &dna[start..end];
            let mut range_emissions = match range.class {
                LabelClass::NonCoding => labeler::noncoding(bases),
                LabelClass::StartCodon => labeler::start_codon(bases),
                LabelClass::StopCodon => labeler::stop_codon(bases, stop_codon_table),
                LabelClass::InternalCodon => labeler::internal_codons(bases, internal_codon_table),
            };
            emissions.append(&mut range_emissions);
        }
        emissions
    }

    pub fn trans_tunings(internal_codon_table: &HashMap<(Base, Base, Base), usize>,
                         stop_codon_table: &HashMap<(Base, Base, Base), usize>)
                         -> Vec<((Label, Label), f64)> {
        let penalty = 5f64;
        let p = 1f64;
        let mut tunings = vec![((Label::NonCoding, Label::StartCodon1), p / penalty),
                               ((Label::StartCodon1, Label::StartCodon2), p),
                               ((Label::StartCodon2, Label::StartCodon3), p)];
        for &id in internal_codon_table.values() {
            tunings.push(((Label::StartCodon3, Label::InternalCodon1(id)), p));
            tunings.push(((Label::InternalCodon1(id), Label::InternalCodon2(id)), p));
            tunings.push(((Label::InternalCodon2(id), Label::InternalCodon3(id)), p));
            for &stop_id in stop_codon_table.values() {
                tunings.push(((Label::InternalCodon3(id), Label::StopCodon1(stop_id)), p));
            }
        }

        for &id in stop_codon_table.values() {
            tunings.push(((Label::StopCodon1(id), Label::StopCodon2(id)), p));
            tunings.push(((Label::StopCodon2(id), Label::StopCodon3(id)), p));
            tunings.push(((Label::StopCodon3(id), Label::NonCoding), p));
        }

        tunings
    }

    pub fn from_labels(start: usize, labels: &Vec<Label>) -> Result<Sequence, String> {
        let mut genes = Vec::new();
        let mut coord_stack = Vec::new();
        for (i, &label) in labels.iter().enumerate() {
            if label == Label::StartCodon1 {
                coord_stack.push(i);
            } else if let Label::StopCodon3(_) = label {
                if let Some(gene_start) = coord_stack.pop() {
                    genes.push((start + gene_start, start + i + 1));
                } else {
                    return Err(format!("found stop codon states outside of gene at {}", i));
                }
            }
        }
        Ok(Sequence {
            start: start,
            end: start + labels.len(),
            genes: genes,
        })
    }

    fn labeled_ranges(&self) -> Vec<LabeledRange> {
        let labeler = |c| {
            move |&r| {
                LabeledRange {
                    class: c,
                    range: r,
                }
            }
        };
        let mut ranges: Vec<LabeledRange> = self.start_codon_ranges()
            .iter()
            .map(labeler(LabelClass::StartCodon))
            .chain(self.stop_codon_ranges()
                .iter()
                .map(labeler(LabelClass::StopCodon))
                .chain(self.internal_codon_ranges()
                    .iter()
                    .map(labeler(LabelClass::InternalCodon))
                    .chain(self.noncoding_ranges()
                        .iter()
                        .map(labeler(LabelClass::NonCoding)))))
            .collect();
        ranges.sort();
        ranges
    }

    fn start_codon_ranges(&self) -> Vec<(usize, usize)> {
        self.genes.iter().cloned().map(|(start, _)| (start, start + 3)).collect()
    }

    fn stop_codon_ranges(&self) -> Vec<(usize, usize)> {
        self.genes.iter().cloned().map(|(_, end)| (end - 3, end)).collect()
    }

    fn internal_codon_ranges(&self) -> Vec<(usize, usize)> {
        self.genes.iter().cloned().map(|(start, end)| (start + 3, end - 3)).collect()
    }

    fn noncoding_ranges(&self) -> Vec<(usize, usize)> {
        let mut noncoding_ranges = Vec::new();
        for i in 0..(self.genes.len()) {
            let (gstart, gend) = self.genes[i];
            if i == 0 {
                noncoding_ranges.push((self.start, gstart));
            } else {
                let (_, last_gend) = self.genes[i - 1];
                noncoding_ranges.push((last_gend, gstart));
            }

            if i == self.genes.len() - 1 {
                noncoding_ranges.push((gend, self.end));
            }
        }
        noncoding_ranges
    }
}

impl FromStr for Sequence {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let terms = s.split_whitespace()
            .map(|s| s.trim_matches(|c| c == ',' || c == '[' || c == ']'))
            .map(|v| v.parse::<usize>().unwrap())
            .collect::<Vec<usize>>();
        let mut pairs = terms.chunks(2).map(|p| (p[0] - 1, p[1]));

        let (seq_start, seq_end) = pairs.next().unwrap();

        let mut genes = Vec::new();
        while let Some((start, end)) = pairs.next() {
            genes.push((start, end))
        }

        Ok(Sequence {
            start: seq_start,
            end: seq_end,
            genes: genes,
        })
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let result = write!(f, "{}-{}:", self.start + 1, self.end);
        let mut gene_results = Vec::new();
        for &(gene_start, gene_end) in self.genes.iter() {
            gene_results.push(write!(f, " [{}, {}]", gene_start + 1, gene_end));
        }
        let composite_result = gene_results.drain(0..).fold(result, |acc, r| acc.and(r));
        composite_result
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn from() {
        assert_eq!("1\t2800\t[190, 255] [337, 2799]".parse::<Sequence>(),
                   Ok(Sequence {
                       start: 0,
                       end: 2800,
                       genes: vec![(189, 255), (336, 2799)],
                   }))
    }

    #[test]
    fn label_dna() {
        let mut internal_codon_table = HashMap::new();
        let mut stop_codon_table = HashMap::new();
        assert_eq!(Sequence {
                           start: 0,
                           end: 14,
                           genes: vec![(1, 13)],
                       }
                       .label_dna(&[Base::A, Base::A, Base::T, Base::G, Base::A, Base::A,
                                    Base::A, Base::C, Base::C, Base::C, Base::T, Base::G,
                                    Base::A, Base::T],
                                  &mut internal_codon_table,
                                  &mut stop_codon_table),
                   vec![(Label::NonCoding, Base::A),
                        (Label::StartCodon1, Base::A),
                        (Label::StartCodon2, Base::T),
                        (Label::StartCodon3, Base::G),
                        (Label::InternalCodon1(0), Base::A),
                        (Label::InternalCodon2(0), Base::A),
                        (Label::InternalCodon3(0), Base::A),
                        (Label::InternalCodon1(1), Base::C),
                        (Label::InternalCodon2(1), Base::C),
                        (Label::InternalCodon3(1), Base::C),
                        (Label::StopCodon1(0), Base::T),
                        (Label::StopCodon2(0), Base::G),
                        (Label::StopCodon3(0), Base::A),
                        (Label::NonCoding, Base::T)])

    }

    #[test]
    fn labeled_ranges() {
        assert_eq!(Sequence {
                           start: 0,
                           end: 1000,
                           genes: vec![(100, 200), (600, 700)],
                       }
                       .labeled_ranges(),
                   vec![LabeledRange {
                            class: LabelClass::NonCoding,
                            range: (0, 100),
                        },
                        LabeledRange {
                            class: LabelClass::StartCodon,
                            range: (100, 103),
                        },
                        LabeledRange {
                            class: LabelClass::InternalCodon,
                            range: (103, 197),
                        },
                        LabeledRange {
                            class: LabelClass::StopCodon,
                            range: (197, 200),
                        },
                        LabeledRange {
                            class: LabelClass::NonCoding,
                            range: (200, 600),
                        },
                        LabeledRange {
                            class: LabelClass::StartCodon,
                            range: (600, 603),
                        },
                        LabeledRange {
                            class: LabelClass::InternalCodon,
                            range: (603, 697),
                        },
                        LabeledRange {
                            class: LabelClass::StopCodon,
                            range: (697, 700),
                        },
                        LabeledRange {
                            class: LabelClass::NonCoding,
                            range: (700, 1000),
                        }]);
    }

    #[test]
    fn noncoding_ranges() {
        assert_eq!(Sequence {
                           start: 0,
                           end: 1000,
                           genes: vec![(100, 200), (600, 999)],
                       }
                       .noncoding_ranges(),
                   vec![(0, 100), (200, 600), (999, 1000)]);
        assert_eq!(Sequence {
                           start: 0,
                           end: 14,
                           genes: vec![(1, 13)],
                       }
                       .noncoding_ranges(),
                   vec![(0, 1), (13, 14)]);
    }

    #[test]
    fn internal_codon_ranges() {
        assert_eq!(Sequence {
                           start: 0,
                           end: 1000,
                           genes: vec![(100, 200), (600, 700)],
                       }
                       .internal_codon_ranges(),
                   vec![(103, 197), (603, 697)])
    }

    #[test]
    fn display() {
        assert_eq!(format!("{}",
                           Sequence {
                               start: 0,
                               end: 2800,
                               genes: vec![(189, 255)],
                           }),
                   "1-2800: [190, 255]")
    }
}
