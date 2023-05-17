use pyo3::prelude::*;
use lj_potential::CalculateLjPotential;
use convert_coordinate::zmat2xyz;
mod lj_potential;
mod convert_coordinate;

#[pymodule]
fn calculate(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(zmat2xyz, m)?)?;
    m.add_class::<CalculateLjPotential>()?;
    Ok(())
}
