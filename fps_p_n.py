import os, sys
import argparse
import subprocess
from clustering import clustering

_HELP_TEXT = """
Usage Example:
    python FPS.py xxx
"""

def get_args():
    parser = argparse.ArgumentParser(description="FPS arguments",
                                     formatter_class=argparse.RawDescriptionHelpFormatter, epilog=_HELP_TEXT)
    parser.add_argument("--msdial_tool_path", help="msdial path")
    parser.add_argument("--deconv_ms2_input_dir", help="input dir path")
    parser.add_argument("--deconv_ms2_output_dir", help="deconvolution output dir path")
    parser.add_argument("--deconv_fps_ms1_input_dir", help="deconvolution fps ms1 input dir path")
    parser.add_argument("--deconv_fps_ms1_output_dir", help="deconvolution fps ms1 input dir path")
    parser.add_argument("--networking_output_dir", help="networking output dir path")
    parser.add_argument("--cluster_output_dir", help='cluster output path')
    parser.add_argument("--msp_path", help="msp file for annotation path")
    parser.add_argument("--annotation_output_dir", help="annotation output dir path")
    args = parser.parse_args()
    return args

def deconv(msdial_path, method_file, input_dir, output_dir, ion_mode):
    res = subprocess.run([msdial_path, 'lcmsdda', '-i', input_dir, '-o',
                          output_dir, '-m', method_file, '-ionMode', ion_mode],
                         stdout=subprocess.PIPE)
    if res.returncode != 0:
        print(res.stdout)
        return False, []
    output_files = []
    # print(res.stdout)
    for result in res.stdout.decode('utf-8').splitlines(True):
        if result.find('msdial_output_file_path') < 0:
            continue
        print(result)
        print(result.split(',')[1].strip())
        output_files.append(result.split(',')[1].strip())
        break
    return True, output_files

def preprocessing(args):
    msdial_path = os.path.join(args.msdial_tool_path, 'MsdialConsoleApp')
    method_file_path = os.path.join(args.deconv_ms2_input_dir, 'msdial.properties')
    ret, deconv_output_files = deconv(msdial_path, method_file_path, args.deconv_ms2_input_dir, args.deconv_ms2_output_dir, 'Positive')
    if not ret:
        print("deconv ms2 failed!")
        return False, None
    deconv_ms2_output_file_path = deconv_output_files[0]
    
    # ms1 deconv
    ms1_method_file_path = os.path.join(args.deconv_ms1_input_dir, 'msdial.properties')
    ret, deconv_ms1_pos_files = deconv(msdial_path, ms1_method_file_path, args.deconv_ms1_input_dir, args.deconv_ms1_output_dir, 'Positive')
    if not ret:
        print("deconv ms1 pos failed!")
        return False, None
    deconv_ms1_pos_output_file_path = deconv_ms1_pos_files[0]
    ret, deconv_ms1_neg_files = deconv(msdial_path, ms1_method_file_path, args.deconv_ms1_input_dir, args.deconv_ms1_output_dir, 'Negative')
    if not ret:
        print("deconv ms1 neg failed!")
        return False, None
    deconv_ms1_neg_output_file_path = deconv_ms1_neg_files[0]
    
    # return True, msdial_output_file_path
    networking_path = os.path.join(args.msdial_tool_path, 'MsdialMolecularNetworkingConsoleApp')
    networking_properties_path = os.path.join(args.deconv_ms2_input_dir, 'networking.properties')
    print("networking_path:" + networking_path)
    print("input:" + deconv_ms2_output_file_path)
    print("output:" + args.networking_output_dir)
    print("network property path:" + networking_properties_path)
    network_res = subprocess.run([networking_path, 'msms', '-i', deconv_ms2_output_file_path,
                                  '-o', args.networking_output_dir, '-p', networking_properties_path],
                                 stdout=subprocess.PIPE)
    # if network_res.returncode != 0:
    #     print(network_res.stdout)
    #     print("networking Error!!!")
    #     return False, None
    netwoking_output_file_path = ''
    for network_result in network_res.stdout.decode('utf-8').splitlines(True):
        print(network_result)
        if network_result.find('candidate_file') < 0:
            continue
        netwoking_output_file_path = network_result.split(',')[1].strip()
        print(netwoking_output_file_path)
        break
    return True, (deconv_ms2_output_file_path, netwoking_output_file_path, deconv_ms1_pos_output_file_path, deconv_ms1_neg_output_file_path)

def fps_pos_neg(args, candidates_file_path, pos_file_path, neg_file_path, deconv_file_path):
    return clustering(args.cluster_output_dir, candidates_file_path, pos_file_path, neg_file_path, deconv_file_path, 2)

def annotation(args, input_file_path):
    annotation_path = os.path.join(args.msdial_tool_path, 'GetAnnotationResult')
    annotation_res = subprocess.run([annotation_path, '--input', input_file_path, '--output', args.annotation_output_dir,
                                     '--library', args.msp_path, '--mz-tolerance', '0.01', '--rt-tolerance', '100'],
                                    stdout=subprocess.PIPE)
    if annotation_res.returncode != 0:
        return False
    return True

def main():
    args = get_args()
    ret, deconv_files = preprocessing(args)
    if not ret or deconv_files is None or len(deconv_files) != 4:
        print("processing error!!!")
        sys.exit(-1)
    print(deconv_files)
    clustering_file_path = fps_pos_neg(args, deconv_files[1], pos_file_path=deconv_files[2],
                neg_file_path=deconv_files[3], deconv_file_path=deconv_files[0])
    # print("clustering_file_path:" + clustering_file_path)
    if not annotation(args, clustering_file_path):
        print("annotation Error!")
    else:
        print("annotation done!!!")

if __name__ == "__main__":
    main()
    