use std::str::FromStr;

#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum Base {
    A,
    C,
    T,
    G,
}

impl Base {
    fn from(c: char) -> Result<Base, String> {
        match c {
            'a' => Ok(Base::A),
            'c' => Ok(Base::C),
            't' => Ok(Base::T),
            'g' => Ok(Base::G),
            x => Err(format!("{} is not a valid dna base", x)),
        }
    }
}

pub fn parse_dna(dna: String) -> Result<Vec<Base>, String> {
    dna.chars().filter(|c| c.is_alphabetic()).map(|c| Base::from(c)).collect()
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_parse() {
        assert_eq!(parse_dna(String::from("ac\nt\ng")),
                   Ok(vec![Base::A, Base::C, Base::T, Base::G]));
        assert_eq!(parse_dna(String::from("ac\nt\nyg")),
                   Err(String::from("y is not a valid dna base")));
    }
}
