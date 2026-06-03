# README — Tutorial Setup Data Custom untuk Ambient Vibration Analysis

Dokumen ini menjelaskan cara menyiapkan data eksperimen sendiri untuk menjalankan script:

```bash
ambient_analysis_user_data_pipeline.py
```

Script ini adalah versi modifikasi dari pipeline MATLAB/Python sebelumnya. Perbedaan utamanya adalah **data sintetis sudah diganti dengan data pengguna** dari file Excel atau CSV. Setelah data terbaca, alur analisis tetap sama:

```text
Data akselerasi eksperimen
        ↓
PSD matrix calculation
        ↓
Singular value spectrum
        ↓
Local modal identification
        ↓
Global mode shape assembly
        ↓
Output frekuensi natural, damping ratio, dan mode shape relatif
```

---

## 1. Instalasi package Python

Jalankan perintah berikut:

```bash
pip install numpy scipy matplotlib pandas openpyxl
```

Jika hanya memakai CSV, `pandas` dan `openpyxl` tidak wajib, tetapi tetap disarankan untuk fleksibilitas pembacaan data.

---

## 2. Struktur file yang disarankan

Letakkan file berikut dalam satu folder kerja:

```text
project_folder/
│
├── ambient_analysis_user_data_pipeline.py
├── ambient_vibration_data_template.xlsx
├── ambient_setup1_example.csv
├── ambient_setup2_example.csv
├── ambient_setup3_example.csv
└── README_ambient_vibration_user_data.md
```

Output akan tersimpan otomatis di folder:

```text
ambient_analysis_output/
```

---

## 3. Format data Excel

Untuk mode Excel, gunakan satu workbook dengan beberapa sheet. Setiap sheet mewakili satu setup sensor.

Contoh:

```text
ambient_vibration_data_template.xlsx
│
├── Setup1
├── Setup2
└── Setup3
```

### Contoh sheet `Setup1`

```csv
time_s,acc_DOF1_mps2,acc_DOF2_mps2,acc_DOF3_mps2,acc_DOF4_mps2
0.000,0.0027,0.0053,0.0079,0.0112
0.005,0.0055,0.0068,0.0111,0.0128
0.010,0.0082,0.0110,0.0135,0.0166
```

### Contoh sheet `Setup2`

```csv
time_s,acc_DOF3_mps2,acc_DOF4_mps2,acc_DOF5_mps2,acc_DOF6_mps2
0.000,0.0041,0.0072,0.0105,0.0128
0.005,0.0064,0.0091,0.0120,0.0145
0.010,0.0078,0.0108,0.0144,0.0170
```

### Contoh sheet `Setup3`

```csv
time_s,acc_DOF2_mps2,acc_DOF4_mps2,acc_DOF5_mps2,acc_DOF6_mps2
0.000,0.0032,0.0075,0.0102,0.0129
0.005,0.0058,0.0094,0.0123,0.0141
0.010,0.0074,0.0115,0.0140,0.0168
```

Kolom pertama **wajib** bernama:

```text
time_s
```

Kolom sinyal akselerasi **wajib** mengikuti format:

```text
acc_DOF1_mps2
acc_DOF2_mps2
acc_DOF3_mps2
...
```

Artinya, data akselerasi harus dalam satuan:

```text
m/s²
```

Jika sensor menghasilkan satuan `g`, konversi dulu:

```text
acc_mps2 = acc_g × 9.80665
```

---

## 4. Format data CSV

Untuk mode CSV, gunakan satu file CSV untuk setiap setup.

Contoh:

```text
ambient_setup1_example.csv
ambient_setup2_example.csv
ambient_setup3_example.csv
```

Isi CSV sama seperti isi sheet Excel. Contoh `ambient_setup1_example.csv`:

```csv
time_s,acc_DOF1_mps2,acc_DOF2_mps2,acc_DOF3_mps2,acc_DOF4_mps2
0.000,0.0027,0.0053,0.0079,0.0112
0.005,0.0055,0.0068,0.0111,0.0128
0.010,0.0082,0.0110,0.0135,0.0166
```

---

## 5. Parameter yang perlu diubah di script

Semua parameter utama ada di bagian:

```python
# =============================================================================
# USER CONFIGURATION
# =============================================================================
```

### 5.1 Pilih mode data

Gunakan Excel:

```python
DATA_MODE = "excel"
```

Atau gunakan CSV:

```python
DATA_MODE = "csv"
```

---

### 5.2 Setting file Excel

Jika memakai Excel:

```python
EXCEL_FILE = "ambient_vibration_data_template.xlsx"
SHEET_NAMES = ["Setup1", "Setup2", "Setup3"]
```

Ubah `EXCEL_FILE` sesuai nama file data Bapak/Ibu.

Contoh:

```python
EXCEL_FILE = "data_getaran_jembatan.xlsx"
SHEET_NAMES = ["Test_1", "Test_2", "Test_3"]
```

Nama sheet harus sama persis dengan nama sheet di Excel.

---

### 5.3 Setting file CSV

Jika memakai CSV:

```python
CSV_FILES = [
    "ambient_setup1_example.csv",
    "ambient_setup2_example.csv",
    "ambient_setup3_example.csv",
]
```

Contoh untuk data sendiri:

```python
CSV_FILES = [
    "uji_getaran_setup1.csv",
    "uji_getaran_setup2.csv",
    "uji_getaran_setup3.csv",
]
```

Urutan file CSV harus sama dengan urutan `SETUP_DOFS`.

---

### 5.4 Setting DOF tiap setup

Contoh default:

```python
SETUP_DOFS = [
    [1, 2, 3, 4],
    [3, 4, 5, 6],
    [2, 4, 5, 6],
]
```

Artinya:

```text
Setup 1 mengukur DOF 1, 2, 3, 4
Setup 2 mengukur DOF 3, 4, 5, 6
Setup 3 mengukur DOF 2, 4, 5, 6
```

Urutan DOF harus sama dengan urutan kolom setelah `time_s`.

Contoh, jika `SETUP_DOFS[0] = [1, 2, 3, 4]`, maka kolom sheet/file Setup 1 harus:

```csv
time_s,acc_DOF1_mps2,acc_DOF2_mps2,acc_DOF3_mps2,acc_DOF4_mps2
```

Jika Setup 1 mengukur DOF 1, 3, 5, dan 7, maka:

```python
SETUP_DOFS = [
    [1, 3, 5, 7],
]
```

Dan kolom datanya harus:

```csv
time_s,acc_DOF1_mps2,acc_DOF3_mps2,acc_DOF5_mps2,acc_DOF7_mps2
```

---

### 5.5 Jumlah DOF global

Default:

```python
N_DOF_GLOBAL = 6
```

Jika struktur punya 8 titik ukur global, ubah menjadi:

```python
N_DOF_GLOBAL = 8
```

Nilai ini harus lebih besar atau sama dengan nomor DOF terbesar dalam `SETUP_DOFS`.

Contoh:

```python
SETUP_DOFS = [
    [1, 2, 3, 4],
    [4, 5, 6, 7],
    [4, 7, 8],
]

N_DOF_GLOBAL = 8
```

---

### 5.6 Reference DOF

Default:

```python
REF_DOF = 4
```

Reference DOF sebaiknya muncul pada semua setup agar mode shape lokal dapat digabungkan lebih stabil.

Contoh setup yang baik:

```python
SETUP_DOFS = [
    [1, 2, 3, 4],
    [4, 5, 6, 7],
    [4, 7, 8, 9],
]
REF_DOF = 4
```

DOF 4 muncul di semua setup.

---

## 6. Parameter frekuensi analisis

### 6.1 Rentang frekuensi analisis PSD dan SV

Default:

```python
F_L = 0.5
F_U = 20.0
MF = 400
```

Artinya analisis dilakukan pada rentang 0.5–20 Hz.

Jika struktur yang diuji diperkirakan memiliki frekuensi natural rendah, misalnya 0.2–10 Hz:

```python
F_L = 0.2
F_U = 10.0
MF = 400
```

Jika struktur kecil atau komponen mesin dengan frekuensi lebih tinggi, misalnya 5–100 Hz:

```python
F_L = 5.0
F_U = 100.0
MF = 800
```

Catatan:

- `F_L` = batas bawah frekuensi analisis.
- `F_U` = batas atas frekuensi analisis.
- `MF` = jumlah interval frekuensi target.
- Semakin besar `MF`, resolusi frekuensi target semakin halus, tetapi data perlu lebih panjang.

---

### 6.2 Band identifikasi mode target

Default:

```python
TARGET_FN_GUESS_HZ = 2.0
TARGET_BAND_LOWER_HZ = 1.6
TARGET_BAND_UPPER_HZ = 2.4
```

Parameter ini dipakai untuk memilih puncak mode yang akan diidentifikasi.

Misalnya dari grafik singular value spectrum terlihat puncak utama di sekitar 3.7 Hz, maka ubah menjadi:

```python
TARGET_FN_GUESS_HZ = 3.7
TARGET_BAND_LOWER_HZ = 3.3
TARGET_BAND_UPPER_HZ = 4.1
```

Jika puncaknya lebar atau damping besar, gunakan band lebih lebar:

```python
TARGET_BAND_LOWER_HZ = 3.0
TARGET_BAND_UPPER_HZ = 4.5
```

Jangan memasukkan dua puncak mode yang terlalu dekat dalam satu band, kecuali memang akan dikembangkan untuk multi-mode identification.

---

## 7. Cara memilih target band dari hasil awal

Langkah praktis:

1. Jalankan script sekali dengan `F_L`, `F_U`, dan `MF` awal.
2. Buka output:

```text
ambient_analysis_output/sv_spectrum_setup1.png
```

3. Cari puncak dominan pada kurva `SV1`.
4. Gunakan puncak tersebut sebagai `TARGET_FN_GUESS_HZ`.
5. Tentukan `TARGET_BAND_LOWER_HZ` dan `TARGET_BAND_UPPER_HZ` agar hanya mencakup satu mode.
6. Jalankan ulang script.

Contoh:

Jika puncak terlihat di 5.8 Hz, ubah:

```python
TARGET_FN_GUESS_HZ = 5.8
TARGET_BAND_LOWER_HZ = 5.3
TARGET_BAND_UPPER_HZ = 6.3
```

---

## 8. Cara menjalankan script

Masuk ke folder kerja, lalu jalankan:

```bash
python ambient_analysis_user_data_pipeline.py
```

Jika berhasil, terminal akan menampilkan informasi seperti:

```text
Loaded user data:
  sampling frequency = 200.0000 Hz
  Setup 1: channels=4, samples=120000, DOFs=[1, 2, 3, 4]
  Setup 2: channels=4, samples=120000, DOFs=[3, 4, 5, 6]
  Setup 3: channels=4, samples=120000, DOFs=[2, 4, 5, 6]

Detected resonance peaks from Setup 1 / SV1:
  2.000 Hz, 5.800 Hz, ...

Local modal identification:
Setup 1: fn=..., zeta=...%
Setup 2: fn=..., zeta=...%
Setup 3: fn=..., zeta=...%
```

---

## 9. Output yang dihasilkan

Output tersimpan di folder:

```text
ambient_analysis_output/
```

Isi output utama:

```text
root_psd_setup1.png
sv_spectrum_setup1.png
global_mode_shape.png
modal_results_summary.csv
global_mode_shape.csv
```

Penjelasan:

| File | Isi |
|---|---|
| `root_psd_setup1.png` | Grafik root PSD untuk Setup 1 |
| `sv_spectrum_setup1.png` | Grafik singular value spectrum untuk mendeteksi puncak mode |
| `global_mode_shape.png` | Grafik mode shape global relatif |
| `modal_results_summary.csv` | Frekuensi natural dan damping ratio tiap setup |
| `global_mode_shape.csv` | Nilai mode shape relatif tiap DOF |

---

## 10. Contoh konfigurasi untuk 1 setup saja

Jika hanya punya satu setup dengan 6 channel:

```python
DATA_MODE = "csv"

CSV_FILES = [
    "data_setup1.csv",
]

SETUP_DOFS = [
    [1, 2, 3, 4, 5, 6],
]

N_DOF_GLOBAL = 6
REF_DOF = 4
```

Format CSV:

```csv
time_s,acc_DOF1_mps2,acc_DOF2_mps2,acc_DOF3_mps2,acc_DOF4_mps2,acc_DOF5_mps2,acc_DOF6_mps2
0.000,...,...,...,...,...,...
0.005,...,...,...,...,...,...
```

Catatan: jika hanya satu setup, hasil mode shape global pada dasarnya sama dengan mode shape lokal pada setup tersebut.

---

## 11. Contoh konfigurasi untuk 8 DOF dan 3 setup

Misalnya struktur punya 8 DOF global dan DOF 4 digunakan sebagai reference:

```python
DATA_MODE = "excel"
EXCEL_FILE = "data_getaran_8dof.xlsx"
SHEET_NAMES = ["Setup1", "Setup2", "Setup3"]

SETUP_DOFS = [
    [1, 2, 3, 4],
    [4, 5, 6, 7],
    [4, 6, 7, 8],
]

N_DOF_GLOBAL = 8
REF_DOF = 4

F_L = 0.5
F_U = 15.0
MF = 500

TARGET_FN_GUESS_HZ = 2.8
TARGET_BAND_LOWER_HZ = 2.4
TARGET_BAND_UPPER_HZ = 3.2
```

---

## 12. Syarat kualitas data eksperimen

Agar hasil lebih stabil:

1. Sampling time harus konstan.
2. Semua setup sebaiknya memiliki sampling frequency yang sama.
3. Satuan semua channel harus sama, disarankan `m/s²`.
4. Arah sensor harus konsisten, misalnya semua arah lateral-X.
5. Durasi data harus cukup panjang.
6. Data tidak boleh memiliki NaN, sel kosong, atau nilai non-numerik.
7. Jika ada trend atau offset besar, script sudah mengurangi nilai rata-rata, tetapi filtering tambahan mungkin tetap diperlukan.
8. Untuk multi-setup, sebaiknya ada reference DOF yang muncul di semua setup.
9. Hindari clipping/saturasi sensor.
10. Hindari data dengan noise dominan jauh lebih besar daripada respons struktur.

Sebagai acuan awal, pipeline asli memakai:

```text
fs = 200 Hz
T = 600 s per setup
```

Jadi jumlah sampel per setup sekitar:

```text
N = 200 × 600 = 120000 sampel
```

Data yang lebih pendek masih bisa diproses, tetapi resolusi frekuensi dan kestabilan identifikasi damping bisa menurun.

---

## 13. Kesalahan umum dan cara memperbaiki

### Error: `Excel file not found`

Penyebab:

- Nama file Excel salah.
- File tidak berada di folder yang sama dengan script.

Solusi:

```python
EXCEL_FILE = "nama_file_yang_benar.xlsx"
```

Atau gunakan path lengkap:

```python
EXCEL_FILE = "D:/Data/Getaran/data_uji.xlsx"
```

---

### Error: `Sheet is missing columns`

Penyebab:

- Nama kolom tidak sesuai format.
- Urutan DOF di `SETUP_DOFS` tidak cocok dengan kolom data.

Contoh, jika:

```python
SETUP_DOFS = [[1, 2, 3, 4]]
```

Maka kolom wajib:

```csv
time_s,acc_DOF1_mps2,acc_DOF2_mps2,acc_DOF3_mps2,acc_DOF4_mps2
```

---

### Warning: `time step is not perfectly uniform`

Penyebab:

- Interval waktu antar sampel tidak konstan.

Solusi:

- Resampling data ke interval waktu tetap.
- Pastikan data logger mengekspor data dengan sampling rate konstan.

---

### Error: `Too few frequency points in target identification band`

Penyebab:

- Band identifikasi terlalu sempit.
- Data terlalu pendek.
- `MF` terlalu kecil atau pengaturan segmentasi tidak cocok.

Solusi:

```python
TARGET_BAND_LOWER_HZ = 1.4
TARGET_BAND_UPPER_HZ = 2.6
```

Atau gunakan data lebih panjang.

---

### Puncak SV tidak terlihat jelas

Kemungkinan penyebab:

- Sensor tidak berada pada lokasi yang sensitif terhadap mode tersebut.
- Data terlalu pendek.
- Excitation ambient terlalu lemah.
- Noise sensor terlalu tinggi.
- Rentang `F_L`–`F_U` tidak sesuai.

Solusi:

- Cek `root_psd_setup1.png` dan `sv_spectrum_setup1.png`.
- Ubah rentang frekuensi analisis.
- Gunakan channel/setup lain yang lebih responsif.
- Tambah durasi pengukuran.

---

## 14. Catatan interpretasi hasil

Untuk data eksperimen, biasanya kita tidak punya `true mode shape`. Karena itu, output seperti `MAC` terhadap mode shape asli dan `RMSE` terhadap mode shape asli tidak tersedia, kecuali Bapak/Ibu punya referensi dari:

- model FEM,
- hasil EMA/impact hammer test,
- hasil pengujian sebelumnya,
- atau hasil simulasi numerik.

Mode shape hasil OMA bersifat relatif. Artinya, besar absolutnya tidak langsung bermakna sebagai deformasi aktual, tetapi pola relatif antar-DOF yang penting.

---

## 15. Checklist sebelum menjalankan data sendiri

Gunakan checklist berikut:

```text
[ ] File Excel/CSV sudah berada di folder yang sama dengan script.
[ ] Kolom pertama bernama time_s.
[ ] Kolom akselerasi mengikuti format acc_DOFx_mps2.
[ ] Satuan akselerasi sudah m/s².
[ ] Sampling time konstan.
[ ] Tidak ada NaN atau sel kosong.
[ ] SETUP_DOFS sesuai dengan kolom data.
[ ] N_DOF_GLOBAL sesuai jumlah DOF total.
[ ] REF_DOF muncul pada semua setup, jika multi-setup.
[ ] F_L dan F_U mencakup frekuensi natural yang dicari.
[ ] TARGET_BAND_LOWER_HZ dan TARGET_BAND_UPPER_HZ hanya mencakup satu mode target.
```

---

## 16. Alur kerja yang disarankan

Untuk penggunaan nyata, gunakan alur berikut:

```text
1. Siapkan data Excel/CSV sesuai format.
2. Jalankan script dengan rentang frekuensi awal yang cukup lebar.
3. Lihat grafik sv_spectrum_setup1.png.
4. Tentukan puncak mode target.
5. Atur TARGET_BAND_LOWER_HZ dan TARGET_BAND_UPPER_HZ.
6. Jalankan ulang script.
7. Cek modal_results_summary.csv.
8. Cek global_mode_shape.csv dan global_mode_shape.png.
9. Jika hasil antar-setup tidak konsisten, periksa sensor, reference DOF, noise, dan target band.
```

---

## 17. Ringkasan parameter penting

| Parameter | Fungsi | Contoh |
|---|---|---|
| `DATA_MODE` | Pilih sumber data | `"excel"` atau `"csv"` |
| `EXCEL_FILE` | Nama file Excel | `"data_getaran.xlsx"` |
| `SHEET_NAMES` | Nama sheet tiap setup | `["Setup1", "Setup2"]` |
| `CSV_FILES` | Nama CSV tiap setup | `["setup1.csv", "setup2.csv"]` |
| `SETUP_DOFS` | DOF yang diukur tiap setup | `[[1,2,3,4],[4,5,6,7]]` |
| `N_DOF_GLOBAL` | Jumlah DOF total | `7` |
| `REF_DOF` | DOF referensi | `4` |
| `F_L` | Batas bawah frekuensi analisis | `0.5` |
| `F_U` | Batas atas frekuensi analisis | `20.0` |
| `MF` | Jumlah interval frekuensi | `400` |
| `TARGET_FN_GUESS_HZ` | Tebakan frekuensi mode target | `2.0` |
| `TARGET_BAND_LOWER_HZ` | Batas bawah band identifikasi | `1.6` |
| `TARGET_BAND_UPPER_HZ` | Batas atas band identifikasi | `2.4` |
| `OUTPUT_DIR` | Folder output | `"ambient_analysis_output"` |

