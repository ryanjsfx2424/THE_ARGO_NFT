## executer.py by Ryan Farber 2022-01-16
"""
The purpose of this script is to modify INPUTS, run the
simulations, and plot.
"""
import os

fn_run  = "nonlinear_advection_periodic_2d.py"
fn_plot = "the_plotter_2d.py"
fn_in   = "INPUTS_NonlinearAdvectionPeriodic2D.py"
fn_in_copy = fn_in + "_copy"
fn_seq = "sequential.py"
fn_mov = "animater.py"

runs = {
         "cfl":2*[0.3,0.29,0.28,0.27,0.26],
        "cmap":["cm.coolwarm", "cm.autumn", "cm.spring", "cm.winter",
                "cm.summer", "cm.Greys", "cm.viridis", "cm.inferno",
                "cm.Reds", "cm.Blues"],
        "type":5*["'explicit'"] + 5*["'implicit'"]
       }

for ii in range(len(runs["cfl"])):
  if runs["type"][ii] == "implicit":
    continue
  # end if

  os.system("cp " + fn_in + " " + fn_in_copy)
  with open(fn_in_copy, "w") as fid_write:
    with open(fn_in, "r") as fid_read:
      for line in fid_read:
        line = line.split()

        if len(line) < 2:
          line = " ".join(line) + "\n"
          fid_write.write(line)
          continue
        # end if

        if   "my_fluid.S" == line[0]:
          line[2] = str(runs["cfl"][ii])
        elif "my_fluid.solver" == line[0]:
          line[2] = str(runs["type"][ii])
        elif "my_fluid.cmap" == line[0]:
          line[2] = str(runs["cmap"][ii])
        # end if/elifs
        line = " ".join(line) + "\n"
        fid_write.write(line)
      # end for
    # end with
  # end with
  os.system("cp " + fn_in_copy + " " + fn_in)
  os.system("python3 " + fn_in)
  os.system("python3 " + fn_run)
  os.system("python3 " + fn_plot)
  os.system("python3 " + fn_seq)
  os.system("python3 " + fn_mov)
# end for ii
