# Mini-Project-Analisis-Pipeline-Sederhana-Sequence-FASTA-FASTQ
Pipeline sederhana untuk menganalisis sequence DNA dari file FASTA/FASTQ

Mini Project Analisis Pipeline Sederhana Sequence FASTA/FASTQ

Mata Kuliah : Struktur Data Bioinformatika

Nama        : Dhafin Fasya Arifin

NIM         : G0401241017

Pipeline ini melakukan:
1. Membaca file FASTA/FASTQ -> disimpan ke List
2. Menghitung frekuensi nukleotida tiap sequence -> pakai Dictionary
3. Mengurutkan sequence berdasarkan GC Content
4. Menampilkan 3 sequence dengan GC Content tertinggi
5. Visualisasi GC Content menggunakan bar chart untuk top 10 sequence tertinggi
6. Menyimpan hasil ke file CSV

```
STRUKTUR FOLDER	

MiniProject/
	├── main.py
	├── data/
	│   └── contoh ecoli.fasta
	├── outputs/
	│   └── grafik_gc_content.png
	│   └── hasil_gc_content.csv
```
