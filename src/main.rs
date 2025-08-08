use std::env;

use tracer::flits::config::FLITSConfig;

fn main() {
    let mut args = env::args().skip(1); // skip program name

    let mode = args
        .next()
        .expect("Expected a mode flag (--flits or --radlite)");
    let value = args.next().expect("Expected a second argument");

    match mode.as_str() {
        "--flits" => {
            let mut flitsc = FLITSConfig::default();
            // use `value` as needed
        }
        "--radlite" => {
            todo!("RadLite functionality not implemented yet");
        }
        _ => {
            eprintln!("Invalid mode: {}", mode);
            std::process::exit(1);
        }
    }

    dbg!(mode, value);
}
