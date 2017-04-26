use std::collections::HashMap;

use corpus::base::Base;

#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum Label {
    NonCoding,
    StartCodon1,
    StartCodon2,
    StartCodon3,
    InternalCodon1(usize),
    InternalCodon2(usize),
    InternalCodon3(usize),
    StopCodon1(usize),
    StopCodon2(usize),
    StopCodon3(usize),
    Invalid,
}

pub fn noncoding(bases: &[Base]) -> Vec<(Label, Base)> {
    bases.iter().cloned().map(|b| (Label::NonCoding, b)).collect()
}

pub fn start_codon(bases: &[Base]) -> Vec<(Label, Base)> {
    vec![(Label::StartCodon1, bases[0]),
         (Label::StartCodon2, bases[1]),
         (Label::StartCodon3, bases[2])]
}

pub fn internal_codons(bases: &[Base],
                       instances: &mut HashMap<(Base, Base, Base), usize>)
                       -> Vec<(Label, Base)> {
    multi_instance_codons(bases, instances, |i, k| match i {
        0 => Label::InternalCodon1(k),
        1 => Label::InternalCodon2(k),
        2 => Label::InternalCodon3(k),
        _ => Label::Invalid,
    })
}

pub fn stop_codon(bases: &[Base],
                  instances: &mut HashMap<(Base, Base, Base), usize>)
                  -> Vec<(Label, Base)> {
    multi_instance_codons(bases, instances, |i, k| match i {
        0 => Label::StopCodon1(k),
        1 => Label::StopCodon2(k),
        2 => Label::StopCodon3(k),
        _ => Label::Invalid,
    })
}

fn multi_instance_codons<F>(bases: &[Base],
                            instances: &mut HashMap<(Base, Base, Base), usize>,
                            f: F)
                            -> Vec<(Label, Base)>
    where F: Fn(usize, usize) -> Label
{
    let mut emissions = Vec::new();
    for codon in bases.chunks(3) {
        let key = (codon[0], codon[1], codon[2]);
        let instance = if let Some(&instance) = instances.get(&key) {
            instance
        } else {
            let instance = instances.len();
            instances.insert(key, instance);
            instance
        };
        for i in 0..3 {
            emissions.push((f(i, instance), codon[i]))
        }
    }
    emissions
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn labeler_internal_codons() {
        let mut unique_checker = HashMap::new();
        unique_checker.insert((Base::A, Base::C, Base::C), 0);
        assert_eq!(internal_codons(&[Base::A, Base::C, Base::C, Base::A, Base::A, Base::A],
                                   &mut unique_checker),
                   vec![(Label::InternalCodon1(0), Base::A),
                        (Label::InternalCodon2(0), Base::C),
                        (Label::InternalCodon3(0), Base::C),
                        (Label::InternalCodon1(1), Base::A),
                        (Label::InternalCodon2(1), Base::A),
                        (Label::InternalCodon3(1), Base::A)]);
    }
}
