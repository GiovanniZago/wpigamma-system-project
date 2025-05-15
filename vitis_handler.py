import os
import h5py
import numpy as np
import subprocess
import config

def get_xilinx_environment():
    command = f'bash -c "source {config.XILINX_SETUP} && env"'
    result = subprocess.run(command, capture_output=True, text=True, shell=True)

    env = {}
    for line in result.stdout.splitlines():
        key, _, value = line.partition("=")
        env[key] = value

    return env

    
def aie_compile_x86(env):
    os.chdir(config.AIE_X86)

    subprocess.run(["aiecompiler", 
                    "--target=x86sim", 
                    "-I", 
                    config.AIE_SRC,
                    config.AIE_SRC + "/graph.cpp", 
                    "-workdir", 
                    config.WORK_X86], 
                    env=env)

def aie_compile_hw(env):
    os.chdir(config.AIE_HW)

    subprocess.run(["aiecompiler", # no debug flag needed here
                    "-target",
                    "hw", 
                    "--platform",
                    f"{config.XILINX_VCK5000_GEN4X8_XDMA}", 
                    "-I", 
                    f"{config.AIE_SRC}", 
                    f"{config.AIE_SRC}/graph.cpp",
                    "-workdir", 
                    f"{config.WORK_HW}", 
                    "--stacksize=16384"], 
                    env=env)
    
def run_x86_simulator(env):
    os.chdir(config.AIE_X86)

    subprocess.run(["x86simulator", 
                    "--dump",
                    f"--pkg-dir={config.WORK_X86}", 
                    f"--input-dir={config.AIE_DATA}",
                    f"--output-dir={config.OUT_SIM_X86}"], 
                    env=env)
    
def run_aiesimulator(env):
    os.chdir(config.AIE_HW)

    subprocess.run(["aiesimulator", 
                    "--dump-vcd=foo", 
                    f"--pkg-dir={config.WORK_HW}",
                    f"--input-dir={config.AIE_DATA}", 
                    f"--output-dir={config.OUT_SIM_AIE}"], 
                    env=env)
    
def hls_compile(env, kernel_name):
    os.chdir(config.HLS)

    subprocess.run(["v++", 
                    "-g", # debug flag
                    "-c", # compile flag
                    "-t", 
                    config.TARGET, 
                    "--platform", 
                    f"{config.XILINX_VCK5000_GEN4X8_XDMA}",
                    "-k",
                    f"{kernel_name}",
                    "-o", 
                    f"{config.HLS}/{kernel_name}/{config.TARGET}/{kernel_name}.{config.TARGET}.xo", 
                    f"{config.HLS}/{kernel_name}/src/{kernel_name}.cpp"], 
                    env=env)

def link_system(env):
    os.chdir(config.LINK)

    aie_dir = "x86" if config.TARGET == "sw_emu" else "hw"

    subprocess.run(["v++", 
                    "-g", # debug flag
                    "-l", # link flag
                    "-t", 
                    config.TARGET, 
                    "--platform", 
                    f"{config.XILINX_VCK5000_GEN4X8_XDMA}", 
                    f"{config.HLS}/mm2s/{config.TARGET}/mm2s.{config.TARGET}.xo", 
                    f"{config.HLS}/s2mm/{config.TARGET}/s2mm.{config.TARGET}.xo", 
                    f"{config.AIE}/{aie_dir}/libadf.a", 
                    "-o", 
                    f"out.{config.TARGET}.xsa", 
                    "--save-temps",
                    "-j", # to speed up computation
                    "4",
                    "--config", 
                    "./system.cfg"], 
                    env=env)

def package_system(env):
    os.chdir(config.PACKAGE)

    aie_dir = "x86" if config.TARGET == "sw_emu" else "hw"

    subprocess.run(["v++",
                    "-g", # debug flag
                    "--package", # package flag
                    "-t",
                    config.TARGET, 
                    "--platform", 
                    f"{config.XILINX_VCK5000_GEN4X8_XDMA}", 
                    f"{config.LINK}/out.{config.TARGET}.xsa", 
                    f"{config.AIE}/{aie_dir}/libadf.a", 
                    "--package.boot_mode=ospi",
                    "-o", 
                    "output.xclbin"], 
                    env=env)

def sw_compile(env):
    os.chdir(config.SW)
    
    if config.TARGET == "sw_emu":
        exec_path = config.SW_EMU

    elif config.TARGET == "hw_emu":
        exec_path = config.HW_EMU

    elif config.TARGET == "hw":
        exec_path = config.HW_RUN

    else:
        print("Please configure target")

    subprocess.run(["g++", 
                    "-std=c++17", 
                    "-g", # debug flag
                    f"-I{config.XRT_INCLUDE}",
                    f"-I{config.HLS_INCLUDE}",
                    "-I./", 
                    f"host_full_dataset.cpp", 
                    "-o", 
                    f"{exec_path}/host.o", 
                    f"-L{config.XRT_LIB}", 
                    f"-L{config.HLS_LIB}",
                    "-lxrt_coreutil", 
                    "-pthread"],
                    env=env)

def run_hw_emu(env):
    os.chdir(config.HW_EMU)

    subprocess.run(["emconfigutil", 
                    "--platform" ,
                    f"{config.XILINX_VCK5000_GEN4X8_XDMA}"], 
                    env=env)
    
    env["XCL_EMULATION_MODE"] = "hw_emu"

    subprocess.run(["./host.o",
                    f"{config.PACKAGE}/output.xclbin"],
                    env=env)
    
def run_hw(env):
    os.chdir(config.HW_RUN)

    subprocess.run(["./host.o",
                    f"{config.PACKAGE}/output.xclbin"],
                    env=env)
    
if __name__ == "__main__":
    env = get_xilinx_environment()

    aie_compile_x86(env)
    run_x86_simulator(env)

    # aie_compile_hw(env)
    # run_aiesimulator(env)

    # hls_compile(env, "mm2s")
    # hls_compile(env, "s2mm")

    # link_system(env)

    # package_system(env)

    # sw_compile(env)

    # run_hw_emu(env)

    # run_hw(env)