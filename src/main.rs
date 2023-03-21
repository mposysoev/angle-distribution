use chemfiles::Frame;
use chemfiles::Trajectory;
use std::env;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufWriter;

fn pbc_distance(atom1: &[f64; 3], atom2: &[f64; 3], box_size: &[f64; 3]) -> f64 {
    let mut delta_x = (atom1[0] - atom2[0]).abs();
    let mut delta_y = (atom1[1] - atom2[1]).abs();
    let mut delta_z = (atom1[2] - atom2[2]).abs();

    if delta_x > box_size[0] / 2.0 {
        delta_x = delta_x - box_size[0];
    }
    if delta_y > box_size[1] / 2.0 {
        delta_y = delta_y - box_size[1];
    }
    if delta_z > box_size[2] / 2.0 {
        delta_z = delta_z - box_size[2];
    }

    return (delta_x * delta_x + delta_y * delta_y + delta_z * delta_z).sqrt();
}

fn cos_alpha(atom1: &[f64; 3], atom2: &[f64; 3], box_size: &[f64; 3]) -> f64 {
    // for any direction
    // let scalar_prod =
    //     (atom1[0] - atom2[0]) * X-component + (atom1[1] - atom2[1]) * Y-component + (atom1[2] - atom2[2]) * Z-component;

    let mut delta_x = atom1[0] - atom2[0];

    if delta_x.abs() > box_size[0] {
        if delta_x > 0.0 {
            delta_x = delta_x - box_size[0];
        } else {
            delta_x = delta_x + box_size[0];
        }
    }

    let scalar_prod = delta_x; // if we take X-axis direction
    let vector_length = pbc_distance(atom1, atom2, box_size);
    return scalar_prod / vector_length;
}

fn angle(atom1: &[f64; 3], atom2: &[f64; 3], box_size: &[f64; 3]) -> Radians {
    return cos_alpha(atom1, atom2, box_size).acos();
}

type Radians = f64;
const CUT_DISTANCE: f64 = 4.0;

fn main() -> std::io::Result<()> {
    let args: Vec<String> = env::args().collect();
    let file_path = &args[1];

    let output_file = File::create(format!("{file_path}-angles.txt"))?;
    let mut writer = BufWriter::new(output_file);

    let mut trajectory = Trajectory::open(file_path, 'r').unwrap();
    let mut frame = Frame::new();

    let file_length = trajectory.nsteps();

    for _i in 0..file_length {
        trajectory.read(&mut frame).unwrap();
        let box_size = frame.cell().lengths();
        let positions = frame.positions();
        for atom in positions {
            for other_atom in positions {
                let dist = pbc_distance(atom, other_atom, &box_size);
                if (dist < CUT_DISTANCE) & !(dist < 0.0001) {
                    let angle_value: Radians = angle(atom, other_atom, &box_size);
                    if !angle_value.is_nan() {
                        writeln!(writer, "{angle_value}")?;
                    }
                }
            }
        }
    }
    Ok(())
}
