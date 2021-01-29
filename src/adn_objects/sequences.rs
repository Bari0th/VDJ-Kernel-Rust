use super::prelude::*;
use std::ops::Deref;

#[derive(Debug)]
pub struct Sequences(Vec<(String, Sequence)>);

impl Sequences {
    pub fn load_from(path: &str) -> Self {
        use std::fs::File;
        use std::io::Read;
        let mut file = File::open(path).expect("Cannot open file");

        let mut buffer = String::new();

        file.read_to_string(&mut buffer).expect("Cannot read file");

        Sequences(
            buffer
                .split('>')
                .skip(1)
                .map(|s| {
                    let mut format_sequence = s.split('\n');
                    (
                        format_sequence.next().unwrap().trim().to_string(),
                        Sequence::new(format_sequence.next().unwrap()),
                    )
                })
                .collect(),
        )
    }
}

impl Deref for Sequences {
    type Target = Vec<(String, Sequence)>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
