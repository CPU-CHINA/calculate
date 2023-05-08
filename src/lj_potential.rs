use std::collections::{BTreeMap, HashMap};
use ordered_float::OrderedFloat;
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyDict, PyList};
use rayon::prelude::*;

#[derive(Eq, PartialEq)]
struct Molecule {
    atoms: BTreeMap<usize, Vec<String>>,
    // bonds: Vec<Vec<usize>>,
    angles: Vec<Vec<usize>>,
    // dihedrals: Vec<Vec<Vec<usize>>>,
    locations: HashMap<String, (OrderedFloat<f64>, OrderedFloat<f64>, OrderedFloat<f64>)>,
}

// 从 PyList 中获取三个元素，并将它们包装成 OrderedFloat<f64> 类型
fn trans_float(num_list: &PyList) -> Result<(OrderedFloat<f64>, OrderedFloat<f64>, OrderedFloat<f64>), PyErr> {
    let first: f64 = num_list.get_item(0).unwrap().extract()?;
    let second: f64 = num_list.get_item(1).unwrap().extract()?;
    let third: f64 = num_list.get_item(2).unwrap().extract()?;

    Ok((
        OrderedFloat(first),
        OrderedFloat(second),
        OrderedFloat(third),
    ))
}

fn create_mole(forest: &PyDict, name: Option<&String>, mole_loc: &PyDict) -> Result<Molecule, PyErr> {
    // 处理 mole_loc
    // ['1MOL', 'C1', '1', [0.147, 1.157, 10.662]] -> HashMap("1": [0.147, 1.157, 10.662])
    let locs: &PyList = mole_loc.get_item(name).unwrap().downcast::<PyList>()?;
    let mut locations: HashMap<String, (OrderedFloat<f64>, OrderedFloat<f64>, OrderedFloat<f64>)> = HashMap::new();
    for mole in locs {
        let mole: &PyList = mole.downcast::<PyList>()?;
        locations.insert(mole.get_item(mole.len() - 2).unwrap().extract()?, trans_float(mole.get_item(mole.len() - 1).unwrap().downcast::<PyList>()?)?);
    }

    let forest_dict = forest.get_item(name).unwrap().downcast::<PyDict>()?;

    let molecule = Molecule {
        atoms: forest_dict.get_item("[atoms]").unwrap().downcast::<PyDict>()?.extract()?,
        // bonds: forest_dict.get_item("[bonds]").unwrap().downcast::<PyList>()?.extract()?,
        angles: forest_dict.get_item("[angles]").unwrap().downcast::<PyList>()?.extract()?,
        // dihedrals: forest_dict.get_item("[dihedrals]").unwrap().downcast::<PyList>()?.extract()?,
        locations,
    };
    Ok(molecule)
}

fn create_forest(moles: &Vec<String>, forest: &PyDict, mole_loc: &PyDict) -> Result<HashMap<String, Molecule>, PyErr> {
    let mut forest_out = HashMap::new();

    for mole in moles {
        let s = mole.to_owned();
        forest_out.insert(s, create_mole(forest, Some(mole), mole_loc)?);
    }

    // No need in calculate LJ potential
    // let complex_dict = forest.get_item("complex").unwrap().downcast::<PyDict>()?;
    // let complex = Molecule {
    //     atoms: complex_dict.get_item("[atoms]").unwrap().downcast::<PyDict>()?.extract()?,
    //     bonds: complex_dict.get_item("[bonds]").unwrap().downcast::<PyList>()?.extract()?,
    //     angles: complex_dict.get_item("[angles]").unwrap().downcast::<PyList>()?.extract()?,
    //     dihedrals: complex_dict.get_item("[dihedrals]").unwrap().downcast::<PyList>()?.extract()?,
    //     locations: BTreeMap::new(),
    // };
    // forest_out.insert(String::from("complex"), complex);

    Ok(forest_out)
}

fn ca_need_atom<'a>(tree_dict: &Vec<Vec<usize>>, atom_set: Vec<&'a usize>, atom: &usize) -> Vec<&'a usize> {
    let mut atom_set = atom_set.clone();
    for angle in tree_dict {
        if angle.contains(&atom) {
            let mut angle = angle.clone();
            angle.retain(|&x| x != *atom);
            atom_set.retain(|&x| !angle.contains(&x));
        }
    }
    atom_set.retain(|&x| x != atom);
    return atom_set;
}

fn get_distance(a: &(OrderedFloat<f64>, OrderedFloat<f64>, OrderedFloat<f64>), b: &(OrderedFloat<f64>, OrderedFloat<f64>, OrderedFloat<f64>)) -> f64 {
    let x = (a.0 - b.0).powi(2);
    let y = (a.1 - b.1).powi(2);
    let z = (a.2 - b.2).powi(2);
    (x + y + z).sqrt()
}

#[pyclass]
pub(crate) struct CalculateLjPotential {
    cc_lib: HashMap<String, Vec<f64>>,
}

#[pymethods]
impl CalculateLjPotential {
    #[new]
    fn new(cc_lib: &PyDict) -> PyResult<Self> {
        let cc_lib: HashMap<String, Vec<f64>> = cc_lib.extract()?;
        Ok(CalculateLjPotential { cc_lib })
    }

    fn calculate(&self, moles: &PyList, forest: &PyDict, mole_loc: &PyDict, lj_cutoff: &PyBool) -> PyResult<f64> {
        let moles = moles.extract::<Vec<String>>()?;
        let forest: HashMap<String, Molecule> = create_forest(&moles, forest, mole_loc)?;
        let lj_cutoff: bool = lj_cutoff.extract()?;

        let e_lj = moles.par_iter().map(|mole| {
            if let Some(molecule) = forest.get(mole) {
                let mut atom_set = molecule.atoms.keys().collect::<Vec<&usize>>();
                let mut e_lj_mole = 0.0;
                for atom1 in molecule.atoms.keys().collect::<Vec<&usize>>() {
                    atom_set = ca_need_atom(&molecule.angles, atom_set, atom1);
                    for atom2 in &atom_set {
                        let a = molecule.locations[&atom1.to_string()];
                        let b = molecule.locations[&atom2.to_string()];
                        let dis = get_distance(&a, &b);
                        if lj_cutoff && dis > 1.4 {
                            continue;
                        } else {
                            let atom_type = vec![&molecule.atoms.get(atom1).unwrap()[1],
                                                 &molecule.atoms.get(atom2).unwrap()[1]];
                            let bond = format!("{}_{}", atom_type[0], atom_type[1]);
                            let c12_c6 = &self.cc_lib.get(&bond).unwrap();
                            e_lj_mole += c12_c6[0] / dis.powi(12) - c12_c6[1] / dis.powi(6);
                        }
                    }
                }
                e_lj_mole
            } else {
                0.0
            }
        }).sum();
        Ok(e_lj)
    }
}
