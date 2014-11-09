#!/usr/bin/python
###########################################################################
# filename      : augment.py
# date written  : 24/08/2013
# written by    : THChew (teonghan@gmail.com)
# description   : read tsv and map locus tag to interproscan id then read
#                 tbl and append matching interproscan ids
# last update   : 24/08/2013    - version 0.1
#                 26/08/2013    - version 0.2 (update db_xref format)
###########################################################################

import sys

if len(sys.argv)!=int(3):
    print "Invalid argument"
    print "Usage: python augment.py tsv_file tbl_file"

else:

    tsv_file=str(sys.argv[1])
    tbl_file=str(sys.argv[2])

    ##tsv_file='LC85_ipr.fa.tsv'
    ##tbl_file='LC85_08192013.tbl'

    import csv
    import collections

    id_map={}
    with open(tsv_file, 'r') as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        for row in csvreader:
            if len(row)>11 and row[11]!='':
                if row[0] not in id_map:
                    id_map[row[0]]='InterPro:%s' % row[11]
                else:
                    id_map[row[0]]='; '.join([id_map[row[0]],('InterPro:%s' % row[11])])

    for i in id_map:
        id_map[i]=list(collections.Counter(id_map[i].split('; ')))
        
    f=open(tbl_file,'r')
    content=f.readlines()
    f.close()

    content=''.join(content)
    content=content.split('>F')

    for i in range(len(content)):
        content[i]=content[i].split('\n')
        for j in range(len(content[i])):
            content[i][j]=content[i][j].split('\t')

    new_content=[]
    for i in range(1,len(content)):
        TMP=[]
        TMP.append('>F%s' % content[i][0][0])
        TMP2=[]
        for j in range(1,len(content[i])):
            if content[i][j][0:3]==['','','']:
                TMP2.append(content[i][j])
            else:
                if len(TMP2)==0:
                    TMP2.append(content[i][j])
                else:
                    TMP.append(TMP2)
                    TMP2=[]
                    TMP2.append(content[i][j])

        new_content.append(TMP)

    del content

    for i in range(len(new_content)):
        TRIGGER=False
        for j in range(1,len(new_content[i])):
            if new_content[i][j][0][-1]=='CDS':
                for k in range(len(new_content[i][j])):
                    if 'locus_tag' in new_content[i][j][k]:
                        locus_tag=new_content[i][j][k][-1]
                        if locus_tag in id_map:
                            TMP=['','','','db_xref']
                            TMP.append(id_map[locus_tag][0]+'\n')
                            for tag in id_map[locus_tag][1:]:
                                TMP.append(2*'\t'+'db_xref\t'+tag+'\n')
                            new_content[i][j].append(TMP)
                            TRIGGER=True

    fout=open('output.tbl','w')
    for i in range(len(new_content)):
        fout.write(new_content[i][0]+'\n')
        for j in new_content[i][1::]:
            for k in j:
                fout.write('\t'.join(k)+'\n')

    fout.close()
