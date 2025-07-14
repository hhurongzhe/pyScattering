import os
import subprocess
import shutil
import sys

files_force = [
    "idahopot.f",
    "av18.f",
    "cdbonn.f",
    "lo450.f",
    "lo500.f",
    "lo550.f",
    "nlo450.f",
    "nlo500.f",
    "nlo550.f",
    "n2lo450.f",
    "n2lo500.f",
    "n2lo550.f",
    "n3lo450new.f",
    "n3lo500new.f",
    "n3lo550new.f",
    "n4lo450.f",
    "n4lo500.f",
    "n4lo550.f",
    "nnlo_opt.f",
    "nnlo_sat.f",
    "n3lo.f",
]


def clean_deps_folder():
    for filename in os.listdir(os.getcwd()):
        if filename.endswith((".o", ".mod", ".dll", ".pyd", ".so", ".dylib", ".pyf")):
            file_path = os.path.join(os.getcwd(), filename)
            os.remove(file_path)


def compile_step1():
    for filename in files_force:
        cmd = ["gfortran", "-c", filename]
        subprocess.run(cmd, check=True)


def compile_step2():
    cmd = "f2py wrapforce.f90 -h wrapforce.pyf -m force_module"
    subprocess.run(cmd, shell=True)


def compile_step3():
    files_step3 = "wrapforce.f90"
    for filename in files_force:
        files_step3 += " " + filename
    cmd = "f2py -c wrapforce.pyf " + files_step3
    print(cmd)
    subprocess.run(cmd, shell=True)


def move_dll_files():
    current_path = os.getcwd()
    libs_folder = os.path.join(current_path, "force_module", ".libs")
    for filename in os.listdir(libs_folder):
        if filename.endswith(".dll"):
            dll_path = os.path.join(libs_folder, filename)
            shutil.copy(dll_path, current_path)

    shutil.rmtree(os.path.join(current_path, "force_module"))


def main():
    clean_deps_folder()
    print("complete clean folders!")
    compile_step1()
    print("complete compiling force package!")
    compile_step2()
    print("complete compiling force module for python!")
    compile_step3()
    print("complete compiling force module for python!")

    if sys.platform == "win32":
        move_dll_files()
        print("complete moving module file to deps folder!")


if __name__ == "__main__":
    main()
