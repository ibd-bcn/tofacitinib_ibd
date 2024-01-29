#!/usr/bin/env bash
#

#Activate environment
#The requirements of this environment are found in requirements.txt
source cpdb/bin/activate

#Cycling_cycling
echo cycling_and_cycling
cd output/cycling_cycling_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo cycling_and_cycling
cd output/cycling_cycling_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo cycling_and_cycling
cd output/cycling_cycling_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo cycling_and_cycling
cd output/cycling_cycling_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#Cycling_epi
echo cycling_and_epi
cd output/cycling_epi_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo cycling_and_epi
cd output/cycling_epi_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo cycling_and_epi
cd output/cycling_epi_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo cycling_and_epi
cd output/cycling_epi_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#Cycling_myeloids
echo cycling_and_myeloids
cd output/cycling_myeloids_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo cycling_and_myeloids
cd output/cycling_myeloids_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo cycling_and_myeloids
cd output/cycling_myeloids_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo cycling_and_myeloids
cd output/cycling_myeloids_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#Cycling_plasmas
echo cycling_and_plasmas
cd output/cycling_plasmas_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo cycling_and_plasmas
cd output/cycling_plasmas_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo cycling_and_plasmas
cd output/cycling_plasmas_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo cycling_and_plasmas
cd output/cycling_plasmas_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#Cycling_stroma
echo cycling_and_stroma
cd output/cycling_stroma_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo cycling_and_stroma
cd output/cycling_stroma_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo cycling_and_stroma
cd output/cycling_stroma_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo cycling_and_stroma
cd output/cycling_stroma_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#Cycling_tcells
echo cycling_and_tcells
cd output/cycling_tcells_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo cycling_and_tcells
cd output/cycling_tcells_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo cycling_and_tcells
cd output/cycling_tcells_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo cycling_and_tcells
cd output/cycling_tcells_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#epi_epi
echo epi_and_epi
cd output/epi_epi_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo epi_and_epi
cd output/epi_epi_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo epi_and_epi
cd output/epi_epi_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo epi_and_epi
cd output/epi_epi_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#epi_plasmas
echo epi_and_plasmas
cd output/epi_plasmas_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo epi_and_plasmas
cd output/epi_plasmas_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo epi_and_plasmas
cd output/epi_plasmas_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo epi_and_plasmas
cd output/epi_plasmas_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#epi_stroma
echo epi_and_stroma
cd output/epi_stroma_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo epi_and_stroma
cd output/epi_stroma_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo epi_and_stroma
cd output/epi_stroma_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo epi_and_stroma
cd output/epi_stroma_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#myeloids_myeloids
echo myeloids_and_myeloids
cd output/myeloids_myeloids_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo myeloids_and_myeloids
cd output/myeloids_myeloids_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo myeloids_and_myeloids
cd output/myeloids_myeloids_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo myeloids_and_myeloids
cd output/myeloids_myeloids_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#myeloids_plasmas
echo myeloids_and_plasmas
cd output/myeloids_plasmas_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo myeloids_and_plasmas
cd output/myeloids_plasmas_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo myeloids_and_plasmas
cd output/myeloids_plasmas_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo myeloids_and_plasmas
cd output/myeloids_plasmas_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#myeloids_stroma
echo myeloids_and_stroma
cd output/myeloids_stroma_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo myeloids_and_stroma
cd output/myeloids_stroma_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo myeloids_and_stroma
cd output/myeloids_stroma_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo myeloids_and_stroma
cd output/myeloids_stroma_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#myeloids_tcells
echo myeloids_and_tcells
cd output/myeloids_tcells_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo myeloids_and_tcells
cd output/myeloids_tcells_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo myeloids_and_tcells
cd output/myeloids_tcells_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo myeloids_and_tcells
cd output/myeloids_tcells_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#plasmas_plasmas
echo plasmas_and_plasmas
cd output/plasmas_plasmas_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo plasmas_and_plasmas
cd output/plasmas_plasmas_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo plasmas_and_plasmas
cd output/plasmas_plasmas_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo plasmas_and_plasmas
cd output/plasmas_plasmas_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#stroma_plasmas
echo stroma_and_plasmas
cd output/stroma_plasmas_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo stroma_and_plasmas
cd output/stroma_plasmas_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo stroma_and_plasmas
cd output/stroma_plasmas_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo stroma_and_plasmas
cd output/stroma_plasmas_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#stroma_stroma
echo stroma_and_stroma
cd output/stroma_stroma_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo stroma_and_stroma
cd output/stroma_stroma_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo stroma_and_stroma
cd output/stroma_stroma_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo stroma_and_stroma
cd output/stroma_stroma_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#tcells_epi
echo tcells_and_epi
cd output/tcells_epi_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo tcells_and_epi
cd output/tcells_epi_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo tcells_and_epi
cd output/tcells_epi_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo tcells_and_epi
cd output/tcells_epi_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#tcells_plasmas
echo tcells_and_plasmas
cd output/tcells_plasmas_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo tcells_and_plasmas
cd output/tcells_plasmas_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo tcells_and_plasmas
cd output/tcells_plasmas_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo tcells_and_plasmas
cd output/tcells_plasmas_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#tcells_stroma
echo tcells_and_stroma
cd output/tcells_stroma_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo tcells_and_stroma
cd output/tcells_stroma_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo tcells_and_stroma
cd output/tcells_stroma_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo tcells_and_stroma
cd output/tcells_stroma_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#tcells_tcells
echo tcells_and_tcells
cd output/tcells_tcells_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo tcells_and_tcells
cd output/tcells_tcells_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo tcells_and_tcells
cd output/tcells_tcells_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo tcells_and_tcells
cd output/tcells_tcells_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#Epi_cycling
echo epi_and_cycling
cd output/epi_cycling_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo epi_and_cycling
cd output/epi_cycling_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo epi_and_cycling
cd output/epi_cycling_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo epi_and_cycling
cd output/epi_cycling_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
#myeloids_epi
echo myeloids_and_epi
cd output/myeloids_epi_group1
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo myeloids_and_epi
cd output/myeloids_epi_group2
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
echo myeloids_and_epi
cd output/myeloids_epi_group3
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10
cd output/myeloids_epi_group4
cellphonedb method statistical_analysis cellphonedb_meta.txt cellphonedb_count.txt --threads=10

