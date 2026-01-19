import sys
from tools.annotate_v2 import VariantAnnotationContainer

if __name__== '__main__':
    if not len(sys.argv) > 1:
        raise NotImplementedError(f"this script takes one or two arguments.\nUsage: annotate_v2.py [sample_name] ([output_path])")
    sample = sys.argv[1]
    if len(sys.argv)>2:
        output_path = sys.argv[2]
    else:
        output_path = f"{sample}.out"
    f = VariantAnnotationContainer(sample, output_path)
    f.print()
