#include <iostream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include "htslib/sam.h"
#include "htslib/thread_pool.h"
using namespace std;


/*
 usage: %prog in.bam thread_num
 g++ -o read_bam_threading -I /path/to/htslib/htslib-1.13 -L /path/to/htslib/htslib-1.13 read_bam_threading.cpp -lhts
 export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/htslib/htslib-1.13
 */
int main(int argc, char* argv[]) {
    int nthread;
    if (argc <= 1) {
        cout << "Please specify a bam file (usage: %prog in.bam thread_num)." << endl;
        return 0;
    } else if (argc <= 2) {
        nthread=1;
        cout << "Set nthread = 1." << endl;
    } else {
        nthread=atoi(argv[2]);
    }
    cout << "Reading " << argv[1] << " with " << nthread << " threads..." << endl;
    
    // open bam
    char* f=argv[1];
    htsFile *in=hts_open(f, "r");
    sam_hdr_t *h=sam_hdr_read(in);
    
    // threading
    htsThreadPool p = {NULL, 0};
    p.pool = hts_tpool_init(nthread);
    hts_set_opt(in, HTS_OPT_THREAD_POOL, &p);
    
    // read bam
    bam1_t *b= bam_init1();
    kstring_t aux={0, 0, NULL};
    while (sam_read1(in, h, b) >= 0) {    // end file = -1
        // details can be found in sam.c bam_read1() and sam.h
        // read
        int      tid      = b->core.tid;                    // target id
        char*     chr     = h->target_name[b->core.tid];    // chr name
        hts_pos_t &start  = b->core.pos;                    // left position, 0-based
        hts_pos_t end     = bam_endpos(b);                  // right position, 0-based
        char*     qname   = bam_get_qname(b);               // read name
        int32_t &l_qseq   = b->core.l_qseq;                 // length of read
        uint16_t &l_qname = b->core.l_qname;                // length of read name
        uint8_t &mapq     = b->core.qual;                   // MAPQ
        uint16_t &flag    = b->core.flag;                   // SAM flag
        uint32_t *cigar   = bam_get_cigar(b);               // cigar, 32-bit
        uint32_t &n_cigar = b->core.n_cigar;                // number of CIGAR operations
        hts_pos_t qlen    = bam_cigar2qlen(n_cigar, cigar); // qlen
        hts_pos_t rlen    = bam_cigar2rlen(n_cigar, cigar); // rlen
        uint16_t &bin     = b->core.bin;                    // bin calculated by bam_reg2bin()
        hts_pos_t &isize  = b->core.isize;                  // insertion size (beween R1 and R2)
        // mate
        char*     mchr    = h->target_name[b->core.mtid];   // mate chr name
        hts_pos_t &mpos   = b->core.mpos;                   // mate left position, 0-based

        // format seq
        uint8_t *tmp_s    = bam_get_seq(b);                 // seq of read, nt16
        char* seq= new char[l_qseq + 1];                    // seq of read, ATGCN
        for (int i=0; i < l_qseq; i++) {
            seq[i]=seq_nt16_str[bam_seqi(tmp_s, i)];        // get nucleotide id and convert into IUPAC id
        }
        seq[l_qseq]='\0';
        
        // format cigar
        string cigarstr;
        for (int k = 0; k < n_cigar; ++k) {
            int leng=bam_cigar_oplen(cigar[k]);
            char opchr=bam_cigar_opchr(cigar[k]);
            cigarstr += to_string(leng) + opchr;
        }
        
        // format query quality
        uint8_t *tmp_q   = bam_get_qual(b);                 // query quality
        char* qqual;
        if (tmp_q[0] == 0xff) {
            qqual= new char[2];
            qqual[0]='*';
            qqual[1]='\0';
        } else {
            qqual= new char[l_qseq + 1];
            for (int i = 0; i < l_qseq; ++i)
                qqual[i]=tmp_q[i]+33;
            qqual[l_qseq]='\0';
        }
        
        // format auxiliary tags
        uint8_t *tmp_a   = bam_get_aux(b);                 // auxiliary data
        uint8_t *tmp_end = b->data + b->l_data;
        vector<char*> auxs;
        while (tmp_a && tmp_end - tmp_a >= 4) {
            aux={0, 0, NULL};
            tmp_a=(uint8_t *)sam_format_aux1(tmp_a, tmp_a[2], tmp_a+3, tmp_end, &aux);
            auxs.push_back(aux.s);
        }
        
        // parse flag
        bool is_paired        = (flag & BAM_FPAIRED) > 0;
        bool is_proper_pair   = (flag & BAM_FPROPER_PAIR) > 0;
        bool is_unmapped      = (flag & BAM_FUNMAP) > 0;
        bool is_mate_unmapped = (flag & BAM_FMUNMAP) > 0;
        bool is_reverse       = (flag & BAM_FREVERSE) > 0;
        bool is_mate_reverse  = (flag & BAM_FMREVERSE) > 0;
        bool is_read1         = (flag & BAM_FREAD1) > 0;
        bool is_read2         = (flag & BAM_FREAD2) > 0;
        bool is_secondary     = (flag & BAM_FSECONDARY) > 0;
        bool is_qc_failed     = (flag & BAM_FQCFAIL) > 0;
        bool is_duplicate     = (flag & BAM_FDUP) > 0;
        bool is_supplementary = (flag & BAM_FSUPPLEMENTARY) > 0;
        
        
        // output
        cout << "tid " << tid << endl;
        cout << "chr " << chr << endl;
        cout << "start " << start << endl;
        cout << "end " << end << endl;
        cout << "qname " << qname << endl;
        cout << "l_qseq " << l_qseq << endl;
        cout << "mapq " << mapq << endl;
        cout << "flag " << flag << endl;
        cout << "cigar " << *cigar << endl;
        cout << "n_cigar " << n_cigar << endl;
        cout << "cigarstr " << cigarstr.c_str() << endl;
        cout << "qlen " << qlen << endl;
        cout << "rlen " << rlen << endl;
        cout << "bin " << bin << endl;
        cout << "mchr " << mchr << endl;
        cout << "mpos " << mpos << endl;
        cout << "isize " << isize << endl;
        cout << "seq " << seq << endl;
        cout << "qqual " << qqual << endl;
        for (char* cp : auxs) {
            cout << "aux " << cp << endl;
        }
        cout << "is_paired " << is_paired << endl;
        cout << "is_proper_pair " << is_proper_pair << endl;
        cout << "is_unmapped " << is_unmapped << endl;
        cout << "is_mate_unmapped " << is_mate_unmapped << endl;
        cout << "is_reverse " << is_reverse << endl;
        cout << "is_mate_reverse " << is_mate_reverse << endl;
        cout << "is_read1 " << is_read1 << endl;
        cout << "is_read2 " << is_read2 << endl;
        cout << "is_secondary " << is_secondary << endl;
        cout << "is_qc_failed " << is_qc_failed << endl;
        cout << "is_duplicate " << is_duplicate << endl;
        cout << "is_supplementary " << is_supplementary << endl;
        cout << endl;
    }
    
    // close bam
    sam_hdr_destroy(h);
    sam_close(in);
    bam_destroy1(b);
    
    // close threads
    hts_tpool_destroy(p.pool);
        
    return 0;
}



/*
 time samtools view test_chrY.bam > /dev/null
 real    0m3.287s
 user    0m3.239s
 sys    0m0.048s

 
 time samtools view -@ 3 test_chrY.bam > /dev/null
 real    0m1.132s
 user    0m3.617s
 sys    0m0.066s
 
 
 time ./read_bam test_chrY.bam > /dev/null
 real    0m3.223s
 user    0m3.127s
 sys    0m0.096s
 
 
 time ./read_bam_threading test_chrY.bam > /dev/null  # 3 threads
 real    0m1.088s
 user    0m3.422s
 sys    0m0.207s
 
 
 time python py_read_bam.py test_chrY.bam > /dev/null
 real    0m3.965s
 user    0m3.941s
 sys    0m0.024s
 
 
 time python py_read_bam.py test_chrY.bam > /dev/null  # 3 threads
 real    0m3.262s
 user    0m4.647s
 sys    0m0.095s

 */
