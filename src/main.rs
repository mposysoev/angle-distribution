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

fn check_pbc(mut delta: f64, box_size: f64) -> f64 {
    if delta.abs() > box_size {
        if delta > 0.0 {
            delta = delta - box_size;
        } else {
            delta = delta + box_size;
        }
    }
    return delta;
}

fn cos_alpha(atom1: &[f64; 3], atom2: &[f64; 3], atom3: &[f64; 3], box_size: &[f64; 3]) -> f64 {
    let len1 = pbc_distance(atom1, atom2, box_size);
    let len2 = pbc_distance(atom1, atom3, box_size);

    let mut delta_x1 = atom2[0] - atom1[0];
    let mut delta_y1 = atom2[1] - atom1[1];
    let mut delta_z1 = atom2[2] - atom1[2];

    let mut delta_x2 = atom3[0] - atom1[0];
    let mut delta_y2 = atom3[1] - atom1[1];
    let mut delta_z2 = atom3[2] - atom1[2];

    delta_x1 = check_pbc(delta_x1, box_size[0]);
    delta_y1 = check_pbc(delta_y1, box_size[1]);
    delta_z1 = check_pbc(delta_z1, box_size[2]);

    delta_x2 = check_pbc(delta_x2, box_size[0]);
    delta_y2 = check_pbc(delta_y2, box_size[1]);
    delta_z2 = check_pbc(delta_z2, box_size[2]);

    let scalar_prod = delta_x1 * delta_x2 + delta_y1 * delta_y2 + delta_z1 * delta_z2;

    return scalar_prod / (len1 * len2);
}

fn angle(atom1: &[f64; 3], atom2: &[f64; 3], atom3: &[f64; 3], box_size: &[f64; 3]) -> Radians {
    return cos_alpha(atom1, atom2, atom3, box_size).acos();
}

type Radians = f64;

fn main() -> std::io::Result<()> {
    let args: Vec<String> = env::args().collect();
    let file_path = &args[1];
    let rcutoff_1_arg = &args[2];
    let rcutoff_2_arg = &args[3];
    let rcutoff_1 = rcutoff_1_arg.parse::<f64>().unwrap(); // 3.8
    let rcutoff_2 = rcutoff_2_arg.parse::<f64>().unwrap(); // 6.2
    println!("Using {rcutoff_1} A cutoff radius for the first shell.");
    println!("Using {rcutoff_2} A cutoff radius for the second shell.");

    let output_file_11 = File::create(format!("{file_path}-angles_11.txt"))?;
    let output_file_12 = File::create(format!("{file_path}-angles_12.txt"))?;
    let output_file_22 = File::create(format!("{file_path}-angles_22.txt"))?;
    let output_file_all = File::create(format!("{file_path}-angles_all.txt"))?;
    let mut writer_11 = BufWriter::new(output_file_11);
    let mut writer_12 = BufWriter::new(output_file_12);
    let mut writer_22 = BufWriter::new(output_file_22);
    let mut writer_all = BufWriter::new(output_file_all);

    let mut trajectory = Trajectory::open(file_path, 'r').unwrap();
    let mut frame = Frame::new();

    let file_length = trajectory.nsteps();

    for _i in 0..file_length {
        trajectory.read(&mut frame).unwrap();
        let box_size = frame.cell().lengths();
        let positions = frame.positions();

        for central_atom in positions {
            for second_atom in positions {
                if central_atom != second_atom {
                    let dist1 = pbc_distance(central_atom, second_atom, &box_size);
                    for third_atom in positions {
                        if third_atom != central_atom && third_atom != second_atom {
                            let dist2 = pbc_distance(central_atom, third_atom, &box_size);
                            if (0.0 < dist1 && dist1 < rcutoff_1) & (0.0 < dist2 && dist2 < rcutoff_1) {
                                let angle_value: Radians =
                                    angle(central_atom, second_atom, third_atom, &box_size);
                                if !angle_value.is_nan() {
                                    writeln!(writer_11, "{angle_value}")?;
                            
                                }
                            }
                            if (0.0 < dist1 && dist1 < rcutoff_1) & (rcutoff_1 < dist2 && dist2 < rcutoff_2) {
                                let angle_value: Radians =
                                    angle(central_atom, second_atom, third_atom, &box_size);
                                if !angle_value.is_nan() {
                                    writeln!(writer_12, "{angle_value}")?;
                            
                                }
                            }
                            if (rcutoff_1 < dist1 && dist1 < rcutoff_2) & (0.0 < dist2 && dist2 < rcutoff_1) {
                                let angle_value: Radians =
                                    angle(central_atom, second_atom, third_atom, &box_size);
                                if !angle_value.is_nan() {
                                    writeln!(writer_12, "{angle_value}")?;
                            
                                }
                            }
                            if (rcutoff_1 < dist1 && dist1 < rcutoff_2) & (rcutoff_1 < dist2 && dist2 < rcutoff_2) {
                                let angle_value: Radians =
                                    angle(central_atom, second_atom, third_atom, &box_size);
                                if !angle_value.is_nan() {
                                    writeln!(writer_22, "{angle_value}")?;
                            
                                }
                            }
                            if (0.0 < dist1 && dist1 < rcutoff_2) & (0.0 < dist2 && dist2 < rcutoff_2) {
                                let angle_value: Radians =
                                    angle(central_atom, second_atom, third_atom, &box_size);
                                if !angle_value.is_nan() {
                                    writeln!(writer_all, "{angle_value}")?;
                            
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(())
}
