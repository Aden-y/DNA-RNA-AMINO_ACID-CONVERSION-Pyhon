import sys
RNA_CODON = {
    'UUU': 'F',     'CUU': 'L',     'AUU': 'I',     'GUU': 'V',
    'UUC': 'F',     'CUC': 'L',     'AUC': 'I',     'GUC': 'V',
    'UUA': 'L',     'CUA': 'L',     'AUA': 'I',     'GUA': 'V',
    'UUG': 'L',     'CUG': 'L',     'AUG': 'M',     'GUG': 'V',
    'UCU': 'S',     'CCU': 'P',     'ACU': 'T',     'GCU': 'A',
    'UCC': 'S',     'CCC': 'P',     'ACC': 'T',     'GCC': 'A',
    'UCA': 'S',     'CCA': 'P',     'ACA': 'T',     'GCA': 'A',
    'UCG': 'S',     'CCG': 'P',     'ACG': 'T',     'GCG': 'A',
    'UAU': 'Y',     'CAU': 'H',     'AAU': 'N',     'GAU': 'D',
    'UAC': 'Y',     'CAC': 'H',     'AAC': 'N',     'GAC': 'D',
    'UAA': 'Stop',  'CAA': 'Q',     'AAA': 'K',     'GAA': 'E',
    'UAG': 'Stop',  'CAG': 'Q',     'AAG': 'K',     'GAG': 'E',
    'UGU': 'C',     'CGU': 'R',     'AGU': 'S',     'GGU': 'G',
    'UGC': 'C',     'CGC': 'R',     'AGC': 'S',     'GGC': 'G',
    'UGA': 'Stop',  'CGA': 'R',     'AGA': 'R',     'GGA': 'G',
    'UGG': 'W',     'CGG': 'R',     'AGG': 'R',     'GGG': 'G'
}
#tc = Codon Translation
#AA = Amino Acid
#co = Codon
#AAS = Amino Acid Sequence (Protein)
#Function to check sequence (split 3) from RNA codon table
#posPro = Possible Protein's Name
#posDNAseq= Possible Protain's DNA sequences
#posRNAseq = possible Protain's RNA sequences
#posRNAseq = possible amino Acid sequences
#s = DNA sequence
def tc(co):#Codon Translation
    AA= None
    if len(co) == 3 and co in RNA_CODON:
        AA = RNA_CODON[co]
    return AA

#Translates The DNA sequence to RNA sequence
def posPro(s):
    ss=[]
    ss=s.replace('T','U')
    results=[]
    AAA = []
    mm=[]
    res=[]
    res1=[]
    l=len(ss)
    for m in range(l):
        AA=tc(ss[m:m+3])
        if AA and AA == 'M':
            mm.append(m)
            found_stop = False
            AAS = ''  # it should only store strings amnino acid sequence
            RA = ''  # store rna sequence of aminoacid sequence above
            DA = ''  # store dna sequence of aminoacid sequence above
            for n in range(m, 1, 3):
                AA = tc(ss[n:n + 3])
                if AA == 'Stop':
                    found_stop=True
                    break
            RA += AAA
            DA += AA
            AAS += AA
            if found_stop:
                results.append(AAS)
                res.append(DA)
                res1.append(RA)
    return (results,res,res1)


if __name__ == "__main__":
    usr=input('Enter the name of the file which contain the DNA sequence with the extension (ex :my.text):: ')
    for i in usr:

        if i == '':
            break
        else:
            parameter= open(usr,'r')
            readingfile=parameter.read()
            possible_a=posPro(readingfile)
            posDNAseq=(possible_a[1])
            posDNAseq1=[]
            for e in posDNAseq:
                if e not in posDNAseq1:
                    posDNAseq1.append(e)
            posRNAseq=(possible_a[2])
            posRNAseq1=[]
            for e in posRNAseq:
                if e not in posRNAseq1:
                    posRNAseq1.append(e)
            posAminoseq=(possible_a[0])
            posAminoseq1=[]
            for e in posAminoseq:
                if e not in posAminoseq1:
                    posAminoseq1.push(e)
            possiblecount=len(posAminoseq1)
            user=input('What Operation would you like to perform?\nDtoR (DNA to RNA) or DtoP (DNA to protein).')
            if user=='DtoR':
                sys.stdout= open('new.txt', 'w')
                for q in range(possiblecount):
                    print ('DNA SEQUENCE: '+' '+posDNAseq1[q]+'\n')
                    print ('RNA SEQUENCE: ' + ' ' + posRNAseq1[q] + '\n\n')
            elif user == 'DtoP':
                sys.stdout = open('new.txt', 'w')
                for q in range(possiblecount):
                    print ('DNA SEQUENCE: ' + ' ' + posDNAseq1[q] + '\n')
                    print ('Amino Acid SEQUENCE: ' + ' ' + posAminoseq1[q] + '\n\n')

            else:
                print ('You have not typed a valid entry. Please run the program agein')