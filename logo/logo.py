# Description -----------------------------------------------------------------
# Create a logo sequence easily
# -----------------------------------------------------------------------------

import argparse
from PIL import Image, ImageDraw, ImageFont
from os.path import realpath, dirname
from re import search

def create_letter(letter, width, height, fnt, color="black", transparency=False):
    "create a letter"
    # font settings
    x_size, y_size = fnt.getsize(letter)
    x_off_set, y_off_set = int(x_size*.02), int(y_size*(.1698 if letter in ["Q", "S"] else 0.2))
    # create image
    img = Image.new('RGBA', (x_size-x_off_set, y_size-y_off_set), 0 if transparency else (255, 255, 255, 255))
    draw = ImageDraw.Draw(img)
    draw.text((-int(float(x_off_set)/2), - y_off_set), letter, color, fnt)
    return img

def resize_logo(logo_pos, width, height):
    "resize the letter"
    try:
        return logo_pos.resize((width, height), Image.ANTIALIAS)
    except:
        return None

def henikoff_pos(full_col, cur_type):
    "Henikoff per position"
    return 1./(len(set(full_col)) * full_col.count(cur_type))

def henikoff_weight(sequences):
    "return a list"
    # reverse seq position
    pos_sequences = ["".join([seq[i] for seq in sequences] for i,_ in enumerate(sequences[0]))]
    henikoff_w = [[] for _ in pos_sequences]
    for pos in pos_sequences:
        for i, el in enumerate(pos):
            henikoff_w[i].append(henikoff_pos(pos, el))
    return [sum(el) for el in henikoff_w]

def get_frequency_pos(sequences, pos, freq=None):
    "compute the frequencies of each aa for a given position"
    results = {}
    for aa in (seq[pos] for seq in sequences):
        if aa in results:
            results[aa] += 1
        else:
            results[aa] = 1
    res_list = [(aa, float(count)/len(sequences)) for aa, count in results.iteritems()]
    res_list.sort(key=lambda el: el[1], reverse=True)
    return res_list

def get_proteus_frequency_pos(sequences, pos, tot=None, freq=None):
    "compute the frequencies of each aa for a given position"
    results = {}
    for aa, count in ((seq[pos], count) for seq, count in sequences):
        if aa in results:
            results[aa] += count
        else:
            results[aa] = count
    res_list = [(aa, float(count)/len(sequences)) for aa, count in results.iteritems()]
    res_list.sort(key=lambda el: el[1], reverse=True)
    return res_list

def get_frequency(sequences, seq_col=3, count_col=5):
    "return a list of frequency per position"
    seq_len = len(sequences[0])
    results = [get_frequency_pos(sequences, pos) for pos in range(seq_len)]
    return results, seq_len

def read_fasta(infile):
    "this function return a list of sequence"
    results, temp_seq = [], ""
    for line in open(infile):
        if line.startswith(">"):
            if temp_seq != "":
                results += [temp_seq]
            temp_seq = ""
        else:
            temp_seq += line.strip()
    if temp_seq != "":
        results += [temp_seq]
    return results

def read_align(infile):
    "this function return a list of sequence from a clustal alignment"
    results, temp_seq = {}, ""
    align_stream = open(infile)
    align_stream.readline()
    for line in align_stream:
        if not line.startswith(" ") and line.strip() != "":
            name, seq = tuple(line.strip().split())
            if name in results:
                results[name] += seq
            else:
                results[name] = seq
    return results.values()

def draw_logo_sequence(logo_img, profile, pos_id, logo_aa, init_top, init_left, pos_width, pos_heigth):
    "draw the logo sequence"
    for pos in pos_id:
        logo_top = init_top
        for aa, freq in profile[pos]:
            pos_logo = resize_logo(logo_aa[aa if aa in logo_aa else "X"], pos_width, int(pos_heigth*freq))
            if pos_logo:
                logo_img.paste(pos_logo, (init_left, logo_top))
                logo_top += pos_logo.size[1]
        init_left += pos_width
    return logo_img

def draw_title(title, img, left_pos, top_pos, x_size, y_size, fnt, transparency):
    "add a title to the left of the logo sequence"
    title = create_letter(title, y_size, x_size, fnt, "black", transparency)
    title = title.rotate(90, expand=1)
    x_t_size, y_t_size = title.size
    x_off_set, y_off_set = (x_t_size - x_size)/2, (y_size - y_t_size)/2
    img.paste(title, (int(x_off_set), int(y_off_set)))
    return img

def draw_aa_logo(color_file, witdh, height, fnt, transparency):
    "return a dictionary that contains color for each type"
    results = {}
    for aa, color in (tuple(line.strip().split()) for line in open(color_file)):
        results[aa] = create_letter(aa, witdh, height, fnt, color, transparency)
    return results

def parse_arguments():
    """Parsing command line
    """
    parser = argparse.ArgumentParser(description="Draw the logo sequence from a seq/fasta/proteus file\nUsage: logo.py <input file> --<seq|fasta|proteus> [options]")
    parser.add_argument('seq_file', default="proteus.dat")
    parser.add_argument('-he', '--heigth', dest='heigth', type=int, default=300)
    parser.add_argument('-wi', '--width', dest='width', type=int, default=1500)
    parser.add_argument('-c', '--color', dest='color_file', default=dirname(__file__)+"/data/color.dat")
    parser.add_argument('-f', '--font', dest='font_file', default=dirname(__file__)+"/data/DejaVuSansMono-Bold.ttf")
    parser.add_argument('-pp', '--position_print', dest='pos_print', nargs="*")
    parser.add_argument('-pi', '--position_id', dest='pos_id', nargs="*", type=int)
    parser.add_argument('-sc', '--seq_col', dest='seq_col', default=3, type=int)
    parser.add_argument('-if', '--id_font', dest='id_font', default=100., type=float)
    parser.add_argument('-o', '--output', dest='output')
    parser.add_argument('-r', '--rotate', action='store_true')
    parser.add_argument('-t', '--title', dest='title')
    parser.add_argument('-tl', '--title_len', dest='title_len', default=.1, type=float)
    parser.add_argument('-tf', '--title_font', dest='title_fnt', default=100., type=float)
    parser.add_argument('--fasta', action='store_true')
    parser.add_argument('--not_show', action='store_false')
    parser.add_argument('--align', action='store_true')
    parser.add_argument('--proteus', action='store_true')
    parser.add_argument('--seq', action='store_true')
    parser.add_argument('--transparency', action='store_true', default=False)
    parser.add_argument('-ra', '--rot_angle', dest='rot_angle', type=int, default=90)
    return parser.parse_args()

def main():
    "draw the logo sequence"
    # Draw args ---------------------------------------------------------------
    args = parse_arguments()
    if args.fasta:
        sequences = read_fasta(args.seq_file)
    elif args.align:
        sequences = read_align(args.seq_file)
    elif args.proteus:
        sequences = [line.split()[args.seq_col] for line in open(args.seq_file) if not line.startswith("#")]
    else:
        sequences = [line.split()[0] for line in open(args.seq_file)]

    assert len(set([len(el) for el in sequences])) == 1, "Sequence(s) of different length!"
    sequence_profile, seq_len = get_frequency(sequences, args.seq_col)
    pos_id = args.pos_id if args.pos_id else [i for i in range(seq_len)]
    pos_print = args.pos_print if args.pos_print else [str(i+1) for i in pos_id]
    assert len(pos_id) == len(pos_print), "Position to print differ from position id!"

    tot_len = len(pos_id)
    script_dir = dirname(__file__)
    fnt = ImageFont.truetype(args.font_file, max([args.heigth, args.width]))

    # Draw sizes --------------------------------------------------------------
    logo_canvas = (args.width, args.heigth)
    logo_img = Image.new('RGBA', logo_canvas, 0 if args.transparency else (255, 255, 255, 255))
    pos_id_heigth = int(float(args.heigth)*0.1)
    title_width = int(args.width * args.title_len) if args.title else 0
    args.width -= title_width
    pos_width, pos_heigth = int(float(args.width - title_width)/tot_len), (args.heigth-pos_id_heigth)
    logo_aa = draw_aa_logo(args.color_file, pos_width, pos_heigth, fnt, args.transparency)
    fnt_pos = ImageFont.truetype(args.font_file, int(float(min(pos_width, pos_id_heigth))/(max(map(len, pos_print))-.5) * (args.id_font/100)))
    fnt_pos_width, fnt_pos_height = fnt_pos.getsize(pos_print[0])

    # Draw logo ---------------------------------------------------------------
    init_top = int(float(args.heigth - ((fnt_pos_width if args.rotate else fnt_pos_height) + pos_heigth))/2)
    init_left = int(float(args.width - pos_width*tot_len)/2)
    init_left += title_width
    logo_im = draw_logo_sequence(logo_img, sequence_profile, pos_id, logo_aa, init_top, init_left, pos_width, pos_heigth)

    # draw position id --------------------------------------------------------
    logo_left = init_left
    for pos in pos_print:
        if args.rotate:
            id_img = create_letter(pos, pos_width, pos_id_heigth, fnt_pos, "black", args.transparency)
            img_y, img_x = id_img.size
            id_img = id_img.rotate(args.rot_angle, expand=1)
        else:
            id_img = create_letter(pos, pos_width, pos_id_heigth, fnt_pos, "black", args.transparency)
            img_x, img_y = id_img.size
        x_off_set = int((pos_width-img_x)*.5)
        y_off_set = int((pos_id_heigth - img_y)*.2)
        logo_img.paste(id_img, (logo_left+x_off_set, pos_heigth+y_off_set+init_top))
        logo_left += pos_width

    if args.title:
        fnt_title = ImageFont.truetype(args.font_file, int(title_width*(args.title_fnt/100.)))
        logo_img = draw_title(args.title, logo_img, init_left - title_width, init_top, title_width, pos_heigth, fnt_title, args.transparency)

    if args.not_show:
        logo_img.show(title=args.title)

    if args.output:
        logo_img.save(args.output, "PNG", dpi=(args.width, args.heigth))

if __name__ == '__main__':
    main()
