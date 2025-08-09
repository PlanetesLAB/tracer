use std::env;
use std::path::Path;

use anyhow::Result;

use tracer::flits::config::FLITSConfig;
use tracer::flits::io::read_fits;

fn main() -> Result<()> {
    let mut args = env::args().skip(1); // skip program name

    let mode = args
        .next()
        .expect("Expected a mode flag (--flits or --radlite)");
    let input_fname = args.next().expect("Expected a second argument");

    match mode.as_str() {
        "--flits" => {
            if args.next().is_some() {
                todo!("Additional arguments for FLITS not supported");
            }
            run_flits(&input_fname)?;
        }
        "--radlite" => {
            todo!("RadLite functionality not implemented yet");
        }
        _ => {
            eprintln!("Invalid mode: {}", mode);
            std::process::exit(1);
        }
    }

    Ok(())
}

fn run_flits<P: AsRef<Path>>(fname: P) -> Result<()> {
    let mut flitsc = FLITSConfig {
        config_path: fname.as_ref().to_path_buf(),
        ..Default::default()
    };
    flitsc.config_dir = flitsc
        .config_path
        .parent()
        .ok_or_else(|| anyhow::anyhow!("Cannot determine config directory"))?
        .to_path_buf();

    flitsc.load_from_path(fname)?;
    read_fits(&mut flitsc)?;
    Ok(())
}
