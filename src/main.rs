use chemfiles::Frame;
use chemfiles::Trajectory;
use rayon::prelude::*;
use std::env;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufWriter;

type Radians = f64;

/// Calculate the distance between two atoms considering periodic boundary conditions.
#[inline]
fn pbc_distance(atom1: &[f64; 3], atom2: &[f64; 3], box_size: &[f64; 3]) -> f64 {
    let mut delta_x = (atom1[0] - atom2[0]).abs();
    let mut delta_y = (atom1[1] - atom2[1]).abs();
    let mut delta_z = (atom1[2] - atom2[2]).abs();

    if delta_x > box_size[0] / 2.0 {
        delta_x -= box_size[0];
    }
    if delta_y > box_size[1] / 2.0 {
        delta_y -= box_size[1];
    }
    if delta_z > box_size[2] / 2.0 {
        delta_z -= box_size[2];
    }

    (delta_x * delta_x + delta_y * delta_y + delta_z * delta_z).sqrt()
}

/// Adjust the delta for periodic boundary conditions.
#[inline]
fn adjust_pbc(mut delta: f64, box_size: f64) -> f64 {
    if delta.abs() > box_size / 2.0 {
        delta -= box_size * delta.signum();
    }
    delta
}

/// Calculate the cosine of the angle between two vectors defined by three atoms.
#[inline]
fn cos_alpha(atom1: &[f64; 3], atom2: &[f64; 3], atom3: &[f64; 3], box_size: &[f64; 3]) -> f64 {
    let len1 = pbc_distance(atom1, atom2, box_size);
    let len2 = pbc_distance(atom1, atom3, box_size);

    let delta_x1 = adjust_pbc(atom2[0] - atom1[0], box_size[0]);
    let delta_y1 = adjust_pbc(atom2[1] - atom1[1], box_size[1]);
    let delta_z1 = adjust_pbc(atom2[2] - atom1[2], box_size[2]);

    let delta_x2 = adjust_pbc(atom3[0] - atom1[0], box_size[0]);
    let delta_y2 = adjust_pbc(atom3[1] - atom1[1], box_size[1]);
    let delta_z2 = adjust_pbc(atom3[2] - atom1[2], box_size[2]);

    let scalar_prod = delta_x1 * delta_x2 + delta_y1 * delta_y2 + delta_z1 * delta_z2;

    scalar_prod / (len1 * len2)
}

/// Calculate the angle between three atoms considering periodic boundary conditions.
#[inline]
fn angle(atom1: &[f64; 3], atom2: &[f64; 3], atom3: &[f64; 3], box_size: &[f64; 3]) -> f64 {
    cos_alpha(atom1, atom2, atom3, box_size).acos()
}

/// Main function to process the trajectory file and calculate angles.
fn main() -> std::io::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: {} <file_path> <rcutoff>", args[0]);
        std::process::exit(1);
    }

    let file_path = &args[1];
    let rcutoff = args[2].parse::<f64>().expect("Invalid cutoff radius");

    println!("Using {} A cutoff radius.", rcutoff);

    let output_file = File::create(format!("{}-angles.txt", file_path))?;
    let mut writer = BufWriter::new(output_file);

    let mut trajectory = Trajectory::open(file_path, 'r').expect("Failed to open trajectory file");
    let mut frame = Frame::new();
    let file_length = trajectory.nsteps();

    let mut all_angles = Vec::new();

    for _ in 0..file_length {
        trajectory.read(&mut frame).expect("Failed to read frame");
        let box_size = frame.cell().lengths();
        let positions = frame.positions().to_vec();

        let frame_angles: Vec<Radians> = calculate_frame_angles(&positions, &box_size, rcutoff);
        all_angles.extend(frame_angles);
    }

    write_angles_to_file(&mut writer, &all_angles)?;

    Ok(())
}

/// Calculate all angles for a single frame.
fn calculate_frame_angles(
    positions: &Vec<[f64; 3]>,
    box_size: &[f64; 3],
    rcutoff: f64,
) -> Vec<Radians> {
    positions
        .par_iter()
        .flat_map(|central_atom| {
            positions
                .par_iter()
                .flat_map(|second_atom| {
                    if central_atom != second_atom {
                        let dist1 = pbc_distance(central_atom, second_atom, box_size);
                        positions
                            .par_iter()
                            .filter_map(move |third_atom| {
                                if third_atom != central_atom && third_atom != second_atom {
                                    let dist2 = pbc_distance(central_atom, third_atom, box_size);
                                    if (0.0 < dist1 && dist1 < rcutoff)
                                        && (0.0 < dist2 && dist2 < rcutoff)
                                    {
                                        let angle_value: Radians =
                                            angle(central_atom, second_atom, third_atom, box_size);
                                        if !angle_value.is_nan() {
                                            return Some(angle_value);
                                        }
                                    }
                                }
                                None
                            })
                            .collect::<Vec<_>>()
                    } else {
                        Vec::new()
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

/// Write calculated angles to the output file.
fn write_angles_to_file(writer: &mut BufWriter<File>, angles: &[Radians]) -> std::io::Result<()> {
    for &angle in angles {
        writeln!(writer, "{}", angle)?;
    }
    Ok(())
}
