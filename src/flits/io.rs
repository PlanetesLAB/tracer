#![allow(unused)]

use anyhow::Result;
use fitsio::FitsFile;

pub fn read_fits(filename: &str) -> Result<()> {
    let mut fptr = FitsFile::open(filename)?;

    let group = 1;
    let first_pix = 1;
    let nullval = -999.0;
    let nullval_d = -999.0;
    // Check dimensions
    let hdu = fptr.primary_hdu()?;
    let mut naxes = [0i32; 3];
    for (i, naxis) in naxes.iter_mut().enumerate() {
        let key_name = format!("NAXIS{}", i + 1);
        *naxis = hdu.read_key::<i32>(&mut fptr, &key_name)?;
    }
    let npixels: i32 = naxes.iter().product();
    let rin: f64 = hdu.read_key(&mut fptr, "Rin")?;
    let rout: f64 = hdu.read_key(&mut fptr, "Rout")?;
    let m_star: f64 = hdu.read_key(&mut fptr, "Mstar")?;
    let r_star: f64 = hdu.read_key(&mut fptr, "Rstar")?;
    let distance: f64 = hdu.read_key(&mut fptr, "distance")?;
    let n_r: i32 = hdu.read_key(&mut fptr, "NXX")?;
    let n_theta: i32 = hdu.read_key(&mut fptr, "NZZ")?;
    let n_spec: i32 = hdu.read_key(&mut fptr, "NSPEC")?;

    todo!()
}
