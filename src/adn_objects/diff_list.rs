use super::prelude::*;
use std::fmt;

pub struct DiffList(Vec<(String, Vec<(String, usize)>)>);

impl DiffList {
    pub fn from_conf(m: &ConfusionMatrix) -> Self {
        let mut diff_list = DiffList(Vec::new());

        let ids = m.get_ids();
        for (idx, id) in ids.iter().enumerate() {
            let diff_line =
                m[idx]
                    .iter()
                    .zip(ids.iter())
                    .fold((Vec::new(), 0), |mut vec, (&cpt, name)| {
                        if cpt > 0 {
                            vec.0.push((name, cpt));
                            vec.1 += cpt;
                        }
                        vec
                    });

            if !diff_line.0.is_empty() {
                let total = diff_line.1;
                let mut v = diff_line
                    .0
                    .into_iter()
                    .map(move |(name, cpt)| (name.clone(), cpt * 100 / total))
                    .collect::<Vec<_>>();
                if m[idx][idx] == 0 {
                    v.push((ids[idx].clone(), 0))
                }
                v.sort_by_key(|(s, cpt)| {
                    if s == &ids[idx] {
                        usize::MIN
                    } else {
                        usize::MAX - cpt
                    }
                });
                diff_list.0.push((ids[idx].clone(), v));
            }
        }

        diff_list
    }
}

impl fmt::Display for DiffList {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for it in self.0.chunks(15) {
            writeln!(
                f,
                "[{}]",
                &it.iter()
                    .map(|(s, _)| format!(", {:>10}", s.strip_prefix("IGH").unwrap()))
                    .collect::<String>()[2..]
            )?;
            writeln!(
                f,
                "{:10?}\n",
                it.iter().map(|(_, v)| v[0].1).collect::<Vec<_>>()
            )?;
        }
        Ok(())
    }
}
