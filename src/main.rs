extern crate rand;
extern crate clap;
extern crate rayon;

use rand::Rng;
use std::{mem, ptr};
use std::time::{Duration, Instant};
use std::num;
use rayon::prelude::*;
use std::sync::{Mutex, Arc, RwLock};
use crossbeam_utils::CachePadded;
use crossbeam_utils::sync::ShardedLock;


const density: f64 = 0.0005;
const mass: f64 = 0.01;
const cutoff: f64 = 0.01;
const min_r: f64 = (cutoff) / 100.0;
const dt: f64 = 0.0005;

const NSTEPS: i32 = 1000;
const SAVEFREQ: i32 = 10;

const n: i32 = 10000;


struct particle_t_accel
{
    ax: f64,
    ay: f64,
    pid: i32,
}


struct particle_t
{
    x: f64,
    y: f64,
    vx: f64,
    vy: f64,
    pid: i32,
}

////
////  keep density constant
////
//fn set_size(n: i32)
//{
//    unsafe {
//        size = (density * n as f64).sqrt();
//    }
//}

//
//  Initialize the particle positions and velocities
//


/* Bin each particle by pos
 *
 *
 * */
fn init_particles(n_t: i32, p_l: &mut [particle_t])
{
//    srand48(time(NULL));

//    let p = p_l.write().unwrap();
    let size: f64 = (density * n as f64).sqrt();

    let mut rng = rand::thread_rng();

    let sx: i32 = ((n as f64).sqrt()).ceil() as i32;
    let sy: i32 = (n + sx - 1) / sx;

    let mut i: usize;
    let mut j: u32;
    let mut k: i32;

    let mut shuffle = vec![0; n as usize];

    for i in 0..n {
        shuffle[i as usize] = i;
    }

    for i in 0..n {

        //
        //  make sure particles are not spatially sorted
        //

        j = (rng.gen::<u32>() % ((n - i) as u32)) as u32; //Generate random number from 0.. n-1
        k = shuffle[j as usize];
        shuffle[j as usize] = shuffle[(n - i - 1) as usize];
        //
        //  distribute particles evenly to ensure proper spacing
        //

        let mut p = &mut p_l[i as usize];

        p.x = size * (1. + (k % sx) as f64) / (1.0 + sx as f64); //gets x
        p.y = size * (1. + (k / sx) as f64) / (1.0 + sy as f64); //gets y


        //
        //  assign random velocities within a bound
        //
        p.vx = (rng.gen::<f64>() * 2.0 - 1.0) as f64;
        p.vy = (rng.gen::<f64>() * 2.0 - 1.0) as f64;
    }
}

fn apply_force(p_acc: &mut particle_t_accel, p: &particle_t, neighbor: &particle_t, dmin: &mut f64, davg: &mut f64, navg: &mut i32)
{
    let dx: f64 = neighbor.x - p.x;
    let dy: f64 = neighbor.y - p.y;

    let mut r2: f64 = dx * dx + dy * dy;
    if r2 > cutoff * cutoff {
        return;
    }
    if r2 != 0.0
    {
        if r2 / (cutoff * cutoff) < *dmin * (*dmin) {
            *dmin = r2.sqrt() / cutoff;
        }
        (*davg) += r2.sqrt() / cutoff;
        (*navg) += 1;
    }

    r2 = r2.max(min_r * min_r);
    let r: f64 = r2.sqrt();

    //
    //  very simple short-range repulsive force
    //

    let coef: f64 = (1.0 - cutoff / r) / r2 / mass;
    p_acc.ax += coef * dx;
    p_acc.ay += coef * dy;
}


//
//  integrate the ODE
//

fn move_particle(p: &mut particle_t, p_acc: &particle_t_accel)
{
    let size: f64 = (density * n as f64).sqrt();
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += p_acc.ax * dt;
    p.vy += p_acc.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    //
    //  bounce from walls
    //
    while p.x < 0.0 || p.x > size {
        if p.x < 0.0 {
            p.x = -1.0 * p.x;
        } else {
            p.x = 2.0 * size - p.x;
        }
        p.vx = -p.vx;
    }
    while p.y < 0.0 || p.y > size {
        if p.y < 0.0 {
            p.y = -1.0 * p.y;
        } else {
            p.y = 2.0 * size - p.y;
        }
        p.vy = -p.vy;
    }
}


fn compute_bin(x: f64, y: f64, size_t: f64, box_num: i32) -> i32 {
    /* Returns bin of particle given its position.
     *
     * */

    return ((y / size_t) * (box_num - 1) as f64 + (x / size_t)).floor() as i32;
}


fn main() {
    let size: f64 = (density * n as f64).sqrt();

    let mut navg = Mutex::new(0);
    let mut davg = Mutex::new(0.0);
    let mut dmin = Mutex::new(1.0);

    //Global stats
    let mut nabsavg: i32 = 0;
    let mut absmin = 1.0;
    let mut absavg = 0.0;

    let mut particles: [particle_t; n as usize] = unsafe { ::std::mem::uninitialized() };
    let mut particles_acc: [particle_t_accel; n as usize] = unsafe { ::std::mem::uninitialized() };

    for i in 0..n {
        particles[i as usize] = particle_t {
            x: 0.0,
            y: 0.0,
            vx: 0.0,
            vy: 0.0,
            pid: i,
        };

        particles_acc[i as usize] = particle_t_accel {
            ax: 0.0,
            ay: 0.0,
            pid: i,
        };
    }

    init_particles(n, &mut particles);

    let sx: i32 = ((n as f64).sqrt()).ceil() as i32;
    let sy: i32 = (n + sx - 1) / sx;

    let old_x: f64 = 0.0;
    let old_y: f64 = 0.0;
    let new_x: f64 = 0.0;
    let new_y: f64 = 0.0;

    let BOX_NUM: i32 = (size / cutoff).floor() as i32;
    let BOX_NUM_t: i32 = BOX_NUM * BOX_NUM;

    let mut k: usize = 0;

    let mut bins = Vec::new();
    let mut bins_read = Vec::new();

    for k in 0..BOX_NUM_t {
//        let mut tmpVec: RwLock<Vec<i32>> = RwLock::new(Vec::new());
//        let mut tmpVec: ShardedLock<Vec<i32>> = ShardedLock::new(Vec::new());
        bins.push(ShardedLock::new(Vec::new()));
        bins_read.push(Vec::new());
    }

    for k in 0..n {
        let bin_index: usize = compute_bin(particles[k as usize].x, particles[k as usize].y, size, BOX_NUM) as usize;
        bins[bin_index].write().unwrap().push(k);
        bins_read[bin_index].push(k);
    }

    let mut step: i32 = 0;
    let mut i: usize = 0;
    let mut j: usize = 0;

    let mut apply_force_time = 0;
    let mut move_time = 0;
    let mut synch_time = 0;

    let mut apply_force_start;
    let mut move_start;
    let mut synch_start;

    let now = Instant::now();

    for step in 0..NSTEPS {

        //
        //  compute forces
        //

        apply_force_start = now.elapsed().as_millis();
        particles_acc.par_iter_mut().for_each(|mut p_acc| {

            //Each thread will have a local copy of the stats
            let mut navg_local = 0;
            let mut davg_local = 0.0;
            let mut dmin_local = 1.0;

            p_acc.ax = 0.0;
            p_acc.ay = 0.0;

            let p = &particles[p_acc.pid as usize];

            let curr_bin: i32 = compute_bin(particles[p_acc.pid as usize].x, particles[p_acc.pid as usize].y, size, BOX_NUM);
            let curr_x: i32 = curr_bin % BOX_NUM;
            let curr_y: i32 = curr_bin / BOX_NUM;

            for j in -1..2 {
                for k in -1..2 {
                    if (curr_x + j >= 0 && curr_x + j < BOX_NUM && curr_y + k >= 0 && curr_y + k < BOX_NUM) {
                        let idx: usize = (curr_x + j + curr_y * BOX_NUM) as usize;
                        bins[idx].read().unwrap().iter()
                                .for_each(|v| {
                                apply_force(p_acc, p, &particles[*v as usize], &mut dmin_local, &mut davg_local, &mut navg_local);
                            }
                            );
                    }
                }
            }


            //Synchronization step
            let mut navg_l = navg.lock().unwrap();
            *navg_l += navg_local;

            let mut davg_l = davg.lock().unwrap();
            *davg_l += davg_local;

            let mut dmin_l = dmin.lock().unwrap();
            if dmin_local < *dmin_l {
                *dmin_l = dmin_local;
            }
        });
        apply_force_time += now.elapsed().as_millis() - apply_force_start;

        //
        //  move particles
        //
        move_start = now.elapsed().as_millis();
        particles.par_iter_mut().for_each(|particle| {
            let old_index: i32 = compute_bin(particle.x, particle.y, size, BOX_NUM);
            move_particle(particle, &particles_acc[particle.pid as usize]);
            let new_index: i32 = compute_bin(particle.x, particle.y, size, BOX_NUM);

            if old_index != new_index {
                {
                    let mut bin = bins[old_index as usize].write().unwrap();
                    bin.retain(|&x| x != particle.pid);
                }
                {
                    let mut bin = bins[new_index as usize].write().unwrap();
                    bin.push(particle.pid)
                }
            }
        });
        move_time += now.elapsed().as_millis() - move_start;

        //
        // Computing statistical data
        //

        synch_start = now.elapsed().as_millis();
        let navg_l = navg.lock().unwrap();
        let davg_l = davg.lock().unwrap();
        let dmin_l = dmin.lock().unwrap();

        if *navg_l > 0 {
            absavg += *davg_l / (*navg_l as f64);
            nabsavg += 1;
        }

        if *(dmin_l) < absmin {
            absmin = *dmin_l;
        }
        synch_time += now.elapsed().as_millis() - synch_start;
    }

    println!("n = {}, simulation time = {} ms, apply_force_time = {}, move_time = {}, synch_time = {}",
             n,
             now.elapsed().as_millis(),
             apply_force_time,
             move_time,
             synch_time);
    
    if nabsavg > 0 {
        absavg /= nabsavg as f64;
    }

    //
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //

    println!(", absmin = {}, absavg = {}", absmin, absavg);
    if absmin < 0.4 {
        println!("The minimum distance is below 0.4 meaning that some particle is not interacting");
    }
    if absavg < 0.8 {
        println!("The average distance is below 0.8 meaning that most particles are not interacting");
    }
}






