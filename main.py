import csv
import os
import matplotlib.pyplot as plt

# folder tempat main.py ini berada
BASE_DIR = os.path.dirname(os.path.abspath(__file__))


#1. Baca file dan simpan ke list
#FASTA
def baca_fasta(filepath):
    list_sequence = []
    seq_id = None
    seq_lines = []

    with open(filepath, "r") as f:
        for baris in f:
            baris = baris.strip()
            if baris == "":
                continue
            if baris.startswith(">"):
                if seq_id is not None:
                    list_sequence.append({
                        "id": seq_id,
                        "sequence": "".join(seq_lines).upper()
                    })
                seq_id = baris[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(baris)

        if seq_id is not None:
            list_sequence.append({
                "id": seq_id,
                "sequence": "".join(seq_lines).upper()
            })

    return list_sequence

#FASTQ
def baca_fastq(filepath):
    list_sequence = []

    with open(filepath, "r") as f:
        baris_list = [b.strip() for b in f if b.strip() != ""]

    for i in range(0, len(baris_list), 4):
        header = baris_list[i]
        sequence = baris_list[i + 1]
        seq_id = header[1:].split()[0]
        list_sequence.append({
            "id": seq_id,
            "sequence": sequence.upper()
        })

    return list_sequence


#2. Hitung frekuensi nukleotida tiap sequence
def hitung_frekuensi_nukleotida(sequence):
    frekuensi = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0}

    for basa in sequence:
        if basa in frekuensi:
            frekuensi[basa] += 1
        else:
            frekuensi["N"] += 1 

    return frekuensi


def hitung_gc_content(frekuensi, panjang_sequence):
    if panjang_sequence == 0:
        return 0.0
    return (frekuensi["G"] + frekuensi["C"]) / panjang_sequence * 100


#Pipeline utama
def main():
    input_file = os.path.join(BASE_DIR, "data", "contoh ecoli.fasta")  #Ganti dengan file FASTA/FASTQ hasil download dari NCBI

    #1. Baca file dan simpan ke list
    os.makedirs(os.path.join(BASE_DIR, "output"), exist_ok=True)
    if input_file.endswith((".fasta", ".fa", ".fna")):
        list_sequence = baca_fasta(input_file)
    elif input_file.endswith((".fastq", ".fq")):
        list_sequence = baca_fastq(input_file)
    else:
        print("Format file tidak didukung. Gunakan file FASTA (.fasta) atau FASTQ (.fastq).")
        return
    print(f"Jumlah sequence terbaca: {len(list_sequence)}\n")

    #2. Hitung frekuensi nukleotida & GC content tiap sequence
    hasil_analisis = []
    for data in list_sequence:
        seq_id = data["id"]
        sequence = data["sequence"]
        panjang = len(sequence)

        frekuensi = hitung_frekuensi_nukleotida(sequence)
        gc_content = hitung_gc_content(frekuensi, panjang)

        hasil_analisis.append({
            "id": seq_id,
            "panjang": panjang,
            "A": frekuensi["A"],
            "T": frekuensi["T"],
            "G": frekuensi["G"],
            "C": frekuensi["C"],
            "GC_content": round(gc_content, 2)
        })

    #3. Urutkan berdasarkan GC content terbesar
    hasil_analisis_sorted = sorted(hasil_analisis, key=lambda x: x["GC_content"], reverse=True)

    #4. Tampilkan 3 sequence dengan GC content tertinggi
    print("=== 3 Sequence dengan GC Content Tertinggi ===")
    for data in hasil_analisis_sorted[:3]:
        print(f"{data['id']:<15} | panjang: {data['panjang']:<5} | GC content: {data['GC_content']}%")

    #5. Visualisasi GC content tiap sequence (bar chart)
    buat_grafik_gc(hasil_analisis_sorted, os.path.join(BASE_DIR, "output", "grafik_gc_content.png"))
    print("\nGrafik disimpan ke output/grafik_gc_content.png")

    #6. Simpan hasil ke CSV
    simpan_ke_csv(hasil_analisis_sorted, os.path.join(BASE_DIR, "output", "hasil_gc_content.csv"))
    print("Hasil disimpan ke output/hasil_gc_content.csv")


def buat_grafik_gc(hasil_analisis_sorted, output_path, top_n=10):
    data_plot = hasil_analisis_sorted[:top_n]
    ids = [d["id"] for d in data_plot]
    gc_values = [d["GC_content"] for d in data_plot]

    plt.figure(figsize=(10, 6))
    warna = ["#2E8B57" if i < 3 else "#4F94CD" for i in range(len(ids))]  # 3 teratas ditandai hijau
    plt.bar(ids, gc_values, color=warna)
    plt.xlabel("Sequence ID")
    plt.ylabel("GC Content (%)")
    plt.title(f"Top {len(ids)} Sequence dengan GC Content Tertinggi")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def simpan_ke_csv(hasil_analisis_sorted, output_path):
    kolom = ["id", "panjang", "A", "T", "G", "C", "GC_content"]

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=kolom)
        writer.writeheader()
        writer.writerows(hasil_analisis_sorted)


if __name__ == "__main__":
    main()
