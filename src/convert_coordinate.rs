use pyo3::prelude::*;
use pyo3::types::{PyDict};
use std::collections::HashMap;


// zmat2cartesian 函数的 Python 接口
#[pyfunction]
pub(crate) fn zmat2xyz(py: Python, zmat: &PyDict) -> PyResult<Py<PyDict>> {
    let mut zmat_hashmap = HashMap::new();
    for (k, v) in zmat.iter() {
        let k = k.extract::<usize>()?;
        let v = v.extract::<[f64; 3]>()?;
        zmat_hashmap.insert(k, v);
    }
    let cartesian = match zmat2cartesian(zmat_hashmap) {
        Ok(cartesian) => cartesian,
        Err(e) => return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(e)),
    };
    let cartesian_dict = PyDict::new(py);
    for (k, v) in cartesian {
        cartesian_dict.set_item(k, v)?;
    }
    Ok(cartesian_dict.into())
}


// 将内坐标转换成笛卡尔坐标
pub(crate) fn zmat2cartesian(zmat: HashMap<usize, [f64; 3]>) -> Result<HashMap<usize, [f64; 3]>, String> {
    let mut cartesian = HashMap::new();
    let mut keys: Vec<&usize> = zmat.keys().collect();
    keys.sort();

    for key in keys {
        let zcoord = zmat.get(key).unwrap();
        let r = zcoord[0];
        let theta = zcoord[1];
        let phi = zcoord[2];

        if *key == 1 {
            cartesian.insert(*key, [0.0, 0.0, 0.0]);
        } else if *key == 2 {
            cartesian.insert(*key, [r, 0.0, 0.0]);
        } else if *key == 3 {
            let x1 = cartesian.get(&(key - 2)).unwrap()[0];
            let x2 = cartesian.get(&(key - 1)).unwrap()[0];
            let z3 = r * theta.cos();
            let y3 = r * theta.sin();
            cartesian.insert(*key, [x1 + (x2 - x1) * z3 / r, y3, 0.0]);
        } else {
            let x2 = cartesian.get(&(key - 2)).unwrap()[0];
            let y2 = cartesian.get(&(key - 2)).unwrap()[1];
            let z2 = cartesian.get(&(key - 2)).unwrap()[2];

            let x3 = cartesian.get(&(key - 1)).unwrap()[0];
            let y3 = cartesian.get(&(key - 1)).unwrap()[1];
            let z3 = cartesian.get(&(key - 1)).unwrap()[2];

            let rx = x3 - x2;
            let ry = y3 - y2;
            let rz = z3 - z2;

            let rxy = (rx.powi(2) + ry.powi(2)).sqrt();

            let cos_phi_32_31_21 = phi.cos();
            let sin_phi_32_31_21 = phi.sin();

            let cos_theta_32_31 = theta.cos();
            let sin_theta_32_31 = theta.sin();

            let x4 =
                x3 + r * (rx * rz * cos_theta_32_31 * cos_phi_32_31_21 / rxy - ry * sin_phi_32_31_21) / r;
            let y4 =
                y3 + r * (ry * rz * cos_theta_32_31 * cos_phi_32_31_21 / rxy + rx * sin_phi_32_31_21) / r;
            let z4 =
                z3 + r * (-rxy * cos_theta_32_31 * cos_phi_32_31_21 + rz * sin_theta_32_31) / r;

            cartesian.insert(*key, [x4, y4, z4]);
        }
    }

    Ok(cartesian)
}