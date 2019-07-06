import pandas as pd

import matplotlib

import tkinter as tk

import allel

import numpy as np

from Bio import SeqIO

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from pyfaidx import Fasta, wrap_sequence
matplotlib.use("TkAgg")



def Split_chromosome():
    print("Enter the FASTA file :")
    fasta = input()
    fasta = '{}.{}'.format(fasta,'fa')
    fa= Fasta(fasta)
    print("Splitting chromosomes..........")

    for seq in fa:
        with open('{}.fa'.format(seq.name), 'w') as out:
            out.write('>{}\n'.format(seq.name))
            for line in wrap_sequence(70, str(seq)):
                out.write(line)
    print("<<<<<<<Splitted>>>>>>")


def chromosome_plotter():

    def onclick(event):
        global Chr, plotted

        i = 0
        for y in ax_y:

            if event.ydata > (y - 0.1) and event.ydata < (y + 0.3):
                Chr = chromosomes[i]
                print("The chromosome is :" + str(chromosomes[i]))

            i += 1

    global chromosomes, vcf

    root=tk.Tk()
    root.title("Variation Visualizer (CHROMOSOME VIEWER)")
    root.state('zoomed')
    df = allel.vcf_to_dataframe(vcf, fields=['CHROM', 'POS', 'REF', 'ALT'], alt_number=1)
    chromosome_vcf = df.CHROM


    fig = Figure(figsize=(10, 20))
    ax = fig.subplots()

    ax.set_ylim(0, 50)
    ax.set_xlim(1, 200)

    canvas = FigureCanvasTkAgg(fig, root)
    fig.set_canvas(canvas=canvas)

    ax_y = []
    ax_labels = []

    p = 0
    m = 48
    v = 47.5

    for chr in chromosomes:
        print(chr)
        chromosome = '{}.{}'.format(chr, 'fa')
        with open(chromosome) as fasta_file:  # Will close handle cleanly
            lengths = []
            for record in SeqIO.parse(fasta_file,
                                      "fasta"):  # (generator)...in this case it is not a multi fasta file
                lengths.append(record.seq)

        seq = str(record.seq)
        print("Initial length :" + str(len(seq)))

        length = int(len(seq) / 2000000)
        print("LEngth is : " + str(length))
        start = 0
        p = 0

        for i in range(start, length):
            ax.broken_barh([(p, 1)], (m, 0.2), facecolors='#1a1a00')
            p += 1

        i = 0
        index = []
        counter = 0
        for chrm in chromosome_vcf:
            if chrm == chr:
                counter += 1
                index.append(i)
            i += 1

        ax.annotate(str(counter) + " variations", (p + 2, m), fontsize=6)
        for j in index:
            pos = df.POS[j]
            pos = int(pos / 2000000)
            ax.broken_barh([(pos, 1)], (v, 1), facecolors='#b30000')

        ax_y.append(m + 0.1)
        ax_labels.append(chr)
        m -= 2
        v -= 2

    plotted = True
    ax.set_yticks(ax_y)
    ax.set_yticklabels(ax_labels, fontsize=12)

    fig.canvas.mpl_connect('button_press_event', onclick)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
    root.mainloop()

def center_window(width=300, height=200):                               # func to set window at centre
                                                                        # get screen width and height
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
                                                                        # calculate position x and y coordinates
    x = (screen_width/2) - (width/2)
    y = (screen_height/2) - (height/2)
    root.geometry('%dx%d+%d+%d' % (width, height, x, y))


def read_fasta():

    global seq, chromosome
    with open(chromosome) as fasta_file:  # Will close handle cleanly
        lengths = []
        for record in SeqIO.parse(fasta_file, "fasta"):  # (generator)...in this case it is not a multi fasta file
            lengths.append(record.seq)

    seq = str(record.seq)
    print("Initial length of the the sequence is :" + str(len(seq)))


def read_vcf():

    global index, df, vcf, Chr
    df = allel.vcf_to_dataframe(vcf, fields=['CHROM', 'POS', 'REF', 'ALT'], alt_number=1)

    chromosome_vcf = df.CHROM
    i = 0
    for chr in chromosome_vcf:
        if chr == Chr:
            index.append(i)
        i += 1
    print(index)
    for j in index:
        print(df.iloc[[j]])


def plot_vcf():

    global position, skipper, checked, df
    for j in index:
        print("POS IS :"+str(df.POS[j]))
        if checked == False:
            a = int(df.POS[j] / skipper)
            position.append(a)

            ax.broken_barh([(a, 4)], (0.56, 1), facecolors='black')
        if checked == True:
            a = int(df.POS[j] /1)
            position.append(a)
            ax.broken_barh([(a, 1)], (0.56, 1), facecolors='black')
    print(position)


def optimiser():
    global start, end, length, skipper, i_length
    length = end - start
    digits = len(str(length))
    print(digits)
    if digits > 3:
        power = digits - 4
        skipper = 2*(10 ** power)
        length = int(length / skipper)
        print("Length is :" + str(length))
        print("Skipper is :" + str(skipper) + "\n\n\n")
        start = 0
        end = length
        i_length = length


def check():
    vcf = 0
    global checked, start, end, skipper, length, position, checked_start, checked_end
    if checked == False:
        if length < 500:

            checked = True
            plot_vcf()
            print("\n\n\nLength is "+str(length))
            print("\n\n\nOTIGINAL SEQ PLOTTED\n\n\n")
            for pos in position:
                if pos > start*skipper and pos < end*skipper:
                    vcf = pos
                    break
            print("The VCF is :" + str(vcf))
            start = vcf - 300
            end = vcf + 300
            checked_start, checked_end = start, end
            length = end - start
            print("Length is "+str(length))
            skipper = 1
            plotter()
            axis_calc()
            canvas.draw()
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)


def plotter():

    global start, end, length, skipper

    i = start
    p = start +1
    m = 10
    counter = 0
    l_counter=0
    length = end - start
    print("Start in plotter :"+ str(start))
    print("end in plotter :" + str(end))

    while counter < length:
        n = seq[i]
        counter += 1
#        print(i)ax.broken_barh([(p, 10)], (0.5, 0.06), facecolors='red')#, hatch='*')
#                                ax.annotate("C", (p + 3, 0.4),fontsize=7)
        if n == 'A' or n == 'a':
            ax.broken_barh([(p, 1)], (0.5, 0.06), facecolors='blue')
            ax.annotate("A", (p + 0.3, 0.4), fontsize=7)
            p += 1

        if n == 'T' or n == 't':
            ax.broken_barh([(p, 1)], (0.5, 0.06), facecolors='yellow')
            ax.annotate("T", (p + 0.3, 0.4), fontsize=7)
            p += 1

        if n == 'G' or n == 'g':
            ax.broken_barh([(p, 1)], (0.5, 0.06), facecolors='green')
            ax.annotate("G", (p + 0.3, 0.4), fontsize=7)
            p += 1

        if n == 'C' or n == 'c':
            ax.broken_barh([(p, 1)], (0.5, 0.06), facecolors='red')
            ax.annotate("C", (p + 0.3, 0.4), fontsize=7)
            p += 1
        if n == 'N':
            ax.broken_barh([(p, 1)], (0.5, 0.06), facecolors='#e6e6e6')
            ax.annotate("N", (p + 0.3, 0.4), fontsize=7)
            p += 1

        i += skipper
    print("\n\n\n"+str(counter)+"   nucleotides plotted.......")


def onclick(event):                                                     # .......function for getting mouse click position

    global start, end, length
    if event.button == 1 :
            print("left click")
            zoomin(event.xdata)
    if event.button == 3 :
            print("right click")
            zoomout(event.xdata)
    if event.button == 2 :
            print("Middle button")
            zoomf()


# ....................................................rescaling axis................................................


def axis_ax1() :
    global Chr
    i=0
    m=0
    print("PLOTTTED>>>>>>>>>")
    ax1.set_ylim(0,1)
    ax1.set_yticks([0.5])
    ax1.set_yticklabels([Chr], fontsize=12)
    ax1.set_xlim(start1,end1)
    ax1.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
    ax1.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
    ax1.set_facecolor('#f2f2f2')
    while i < length:
        ax1.broken_barh([(m, 1)], (0.5, 0.2), facecolors='#4d4d4d')
        m += 1
        i += 1


def axis_calc() :
    global Chr
    ax.set_ylim(0, 1.2)
    ax.set_yticks([0.41])

    ax.set_xlim(start,end)
    if checked == False:
        ax.set_xlabel('Hypothetical positions of '+Chr, fontsize=12)
        ax.set_yticklabels(['Reference seq.'], fontsize=12)
    else:
        ax.set_xlabel('Actual positions of '+Chr, fontsize=12)
        ax.set_yticklabels(['Actual seq.'], fontsize=12)
    ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)


prev_i,prev_j = 0,0


def miniature_view():
    m_i=int(start)
    m_j=int(end)

    global prev_i,prev_j, skipper
    while prev_i < prev_j :
        ax1.broken_barh([(prev_i, 1)], (0.5, 0.2), facecolors='#595959')
        prev_i+=1
    prev_i=m_i
    prev_j=m_j
    print ("miniture start is" + str(m_i) )
    print("\nminiture end is" + str(m_j))
    while m_i < m_j :
        ax1.broken_barh([(m_i, 1)], (0.5, 0.2), facecolors='white')
        m_i+=1
    for pos in index:
        a = int(df.POS[pos] / skipper)
        ax1.broken_barh([(a, 4)], (0, 1), facecolors='red')
    canvas1.draw()
    canvas1.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)


def zoomf() :
    global start,end,prev_positiontozoom,positiontozoom
    start=0
    end=i_length
    ax.broken_barh([(prev_positiontozoom, 10)], (0.3, 0.02), facecolors='white')
    axis_calc()
    miniature_view()
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    prev_positiontozoom = positiontozoom


def zoomin(positiontozoom):
    global end, start, length, prev_positiontozoom
    length=end - start
    start = int(positiontozoom - 0.3*length)

    end = int(positiontozoom + 0.3*length)

    if checked == False:
        miniature_view()
        if start < 0:
            start = 0
        if end > i_length:
            end = i_length
        check()
        ax.broken_barh([(prev_positiontozoom, 2)], (0.3, 0.02), facecolors='white')
        prev_positiontozoom = positiontozoom
    if checked == True :
        ax.broken_barh([(prev_positiontozoom, 2)], (0.3, 0.02), facecolors='white')
        prev_positiontozoom = positiontozoom
    axis_calc()
    ax.broken_barh([(positiontozoom, 2)], (0.3, 0.02), facecolors='black')
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)


def zoomout(positiontozoom):
    global end, start, length, prev_positiontozoom
    length = end - start
    start = int(start - 0.3 * length)

    end = int(end + 0.3 * length)

    if checked == False:
        miniature_view()
        if start < 0:
            start = 0
        if end > i_length:
            end = i_length
        ax.broken_barh([(prev_positiontozoom, 2)], (0.3, 0.02), facecolors='white')
        prev_positiontozoom = positiontozoom
    if checked == True:
        ax.broken_barh([(prev_positiontozoom, 2)], (0.3, 0.02), facecolors='white')
        prev_positiontozoom = positiontozoom
    axis_calc()
    ax.broken_barh([(positiontozoom, 2)], (0.3, 0.02), facecolors='black')
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)


def update_annot(x,y,counter):

    annot.xy = (x,y)
    text=''
    text = "Chromosome : "+df.CHROM[counter]+"\nPosition  : "+str(df.POS[counter])+"\nREF         : "+df.REF[counter]+"\nALT         : "+df.ALT[counter]
    annot.set_text(text)
    annot.get_bbox_patch().set_facecolor("red")
    annot.get_bbox_patch().set_alpha(0.4)


def hover(event):
    x, y = event.xdata, event.ydata

    global position, index

    if y.__class__ == np.float64 :

        if y > 0.56 and y < 0.8 :
            vis = annot.get_visible()
            counter = index[0]
            x = int(x)
            for pos in position :
                if x >= int(pos) and x < int(pos) + 4 :


                    if checked == True:
                        print("Hovering over variation no. : " + str(pos))
                        for i in index:
                            if df.POS[i] == pos:
                                counter = i
                        update_annot(x, y, counter)
                    else :
                        print("Hovering over variation no. : " + str(2 * pos))

                        update_annot(x, y, counter)

                    annot.set_visible(True)
                    fig.canvas.draw_idle()
                else :
                    if vis :
                        annot.set_visible(False)
                        fig.canvas.draw_idle()
                counter = counter + 1


def getit(event):                        # ............function to detect if slider is increasing or decresing
    zoom1 = w.get()
    global zoom2,zoom,start,end

    scale = 1.0
    if zoom1 == 0:
        zoom()
    if zoom2 < zoom1:
        print("increasing\n\n\n")
        zoom2 = zoom1
        mid=start+((end-start)/2)
        zoomin(mid)

    else:
        print("decreasing\n\n\n")
        zoom2 = zoom1
        mid = start + ((end - start) / 2)
        zoomout(mid)


    print(zoom2)


def go_left():
    global start, end,length, checked_start, min_count

    min_count += 1

    length = end - start

    if checked:

        start = int(start - (0.1*length))
        end = int(end - (0.1*length))

        if start < checked_start:
            print("\n\nPLOTTED  AGAIN......\n\n")
            plotter()
    else:

        start = int(start - (0.01 * length))
        end = int(end - (0.01 * length))
        if start < 0:
            print("OUT OF BOUNDARY")
            start = 0
        if min_count ==5:
            min_count = 0
            miniature_view()

    axis_calc()
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)


def go_right():
    global start, end, length, checked_end, min_count
    min_count += 1
    length = end - start

    if checked:
        start = int(start + (0.1 * length))
        end = int(end + (0.1 * length))

        if end > checked_end:
            print("\n\nPLOTTED  AGAIN......\n\n")
            plotter()
    else:
        start = start + int(0.01 * length)
        end = end + int(0.01 * length)
        if end > (start + length):
            end = start + length
            print("OUT OF BOUNDARY")
        if min_count ==5:
            min_count = 0
            miniature_view()

    axis_calc()
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)


# .............................setting tkinter window..........................

print("Do you want to split the Genome file?")
ans = input()

if ans == 'y' or ans == 'Y' or ans == 'yes' or ans == 'Yes':
    Split_chromosome()

print("Enter the VCF file :")
vcf = input()
vcf = '{}.{}'.format(vcf,'vcf')
print(vcf)

chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                   'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
                   'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
Chr=''
chromosome_plotter()
print("Now Chr is :"+str(Chr))
chromosome = '{}.{}'.format(Chr,'fa')
print(chromosome)


root = tk.Tk()
root.state('zoomed')
root.title("Variation Visualizer")
prev_positiontozoom = 0

zoom1, zoom2, zoom = 1, 1, 0

w = tk.Scale(root, length=300, width=20, sliderlength=10, tickinterval=2, from_=0, to=10,resolution=2,
             orient=tk.HORIZONTAL,command = getit)

w.pack(side=tk.BOTTOM)


# .............................setting up matplotlib figure and plot................................


fig = Figure(figsize=(10,5))

fig1 = Figure(figsize=(10,1))

ax = fig.subplots()
ax1 = fig1.subplots()

canvas=FigureCanvasTkAgg(fig,root)
canvas1=FigureCanvasTkAgg(fig1,root)
fig.set_canvas(canvas=canvas)
fig1.set_canvas(canvas=canvas1)

button1 = tk.Button(root, text = "<-- LEFT" , command = go_left)
button1.configure(width = 10, activebackground = "#bfbfbf", relief = tk.FLAT)

button1.pack(side=tk.LEFT)
button2 = tk.Button(root, text = "RIGHT-->" ,command = go_right)
button2.configure(width = 10, activebackground = "#bfbfbf", relief = tk.FLAT)
button2.pack(side=tk.RIGHT)


# ...................................Hover .........................................


annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->"))
annot.set_visible(False)

# ...............................fasta, vcf and variables...................................

seq = ''
index = []
position = []
df = pd.DataFrame()
checked = False
checked_start, checked_end =0, 0
min_count = 0


read_fasta()

i_length = len(seq)

length, start, end = 0, 0, i_length          # ........Global variables declared

skipper = 0



optimiser()

plotter()

end1 = i_length + 500
start1 = -500

axis_calc()
axis_ax1()

read_vcf()
plot_vcf()

miniature_view()

fig.canvas.mpl_connect('button_press_event', onclick)                   # ......event (  click  )
fig.canvas.mpl_connect("motion_notify_event", hover)                    # ......event ( hover   )

# ...............plot................

axis_calc()
axis_ax1()
canvas.draw()
canvas1.draw()
canvas1.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

root.mainloop()