use pyo3::prelude::*;
use lj_potential::CalculateLjPotential;

mod lj_potential;

#[pymodule]
fn calculate(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<CalculateLjPotential>()?;
    Ok(())
}
