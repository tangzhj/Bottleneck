# scATAC-seq analyses

## Input files
 Documenting in ./Input data<br>
 **annovar_result**:The mitochondrial variants were annotated with ANNOVAR.The annotated variants include loop (loop mutation), tRNA (tRNA mutation), rRNA (rRNA mutation), coding (coding-region mutation), NS (non-synonymous mutation), and SY (synonymous mutation) according to variants locations
 <br>**matrix**: Matrix of variant in each cell type(common consist of file suffix with '_matrix_count.txt', '_matrix_variation.txt', '_Site_infor.txt')
 <br>**cell_mutation_count.txt**: Mutation count in each cell(annotation of cell type).
 <br>**cell_num_with_variation.txt**: Mutation count in each cell cellypte.
 <br>**mono_ratio.txt**: Variant informations shared in progenitor cell and myeloid cell.(name:variant;a/b:cell num with VAF<0.5/VAF>0.5 in progenitor cell;<br>c/d:cell num with VAF<0.5/VAF>0.5 in myeloid cell;)
 <br>**tbnk_ratio.txt**: Variant informations shared in progenitor cell and lymphoid cell.(name:variant;a/b:cell num with VAF<0.5/VAF>0.5 in progenitor cell;<br>c/d:cell num with VAF<0.5/VAF>0.5 in lymphoid cell;)
 <br>**proj_res1_annotation_result.txt**: mtDNA reads in each annotated cell.
 <br>**cell_mutation_count.txt**:
## Analysis content
 1.The number of somatic mtDNA mutations per cell. Figure. 1c.<br>
 2.Allele frequency spectrum of somatic mtDNA mutations. Figure. 1d.<br>
 3.relative mtDNA copies in mtscATAC-seq data. Figure. 1e.<br>
 4.Signals of purifying selection at specific mtDNA genomic site. Figure.4<br>
 5.%cell with variations. Figure. S1b.<br>
 6.The allele frequency spectrum of somatic mtDNA mutations for different types. Figure. S4.<br>
 
