for coupling in L12 L13 L21 L22 L23 L31 L32 L33 R11 R12 R13 R21 R22 R23 R31 R32 R33
do
    scp -r 10.2.12.77:/home/atirek/Documents/MG5_aMC_v3_2_0/cross_section/"$coupling"/interference/results.txt "$coupling"/interference.txt
    scp -r 10.2.12.77:/home/atirek/Documents/MG5_aMC_v3_2_0/cross_section/"$coupling"/pair/results.txt "$coupling"/pair.txt
    scp -r 10.2.12.77:/home/atirek/Documents/MG5_aMC_v3_2_0/cross_section/"$coupling"/pureqcd/results.txt "$coupling"/pureqcd.txt
    scp -r 10.2.12.77:/home/atirek/Documents/MG5_aMC_v3_2_0/cross_section/"$coupling"/single/results.txt "$coupling"/single.txt
    scp -r 10.2.12.77:/home/atirek/Documents/MG5_aMC_v3_2_0/cross_section/"$coupling"/t_channel/results.txt "$coupling"/t_channel.txt
done