# MetaboDesk Karşılaştırmalı Doğrulama Planı

## 🎯 Amaç
MetaboDesk'in ürettiği sonuçların, geleneksel komut satırı araçlarıyla (COBRApy, memote CLI) birebir aynı olduğunu kanıtlamak.

## 🧬 Test Modeli
**E. coli Core Model** (`e_coli_core.xml`)
- 72 metabolit, 95 reaksiyon, 137 gen
- Biyoinformatik alanında en yaygın kullanılan benchmark model
- COBRApy ile birlikte gelir (`cobra.io.load_json_model("textbook")` veya doğrudan SBML)

---

## 📋 Test Adımları (Toplam 10 Test)

### FAZA 1: Model Yükleme & Bilgi Doğrulama
| # | Test | Geleneksel Araç (Python/CLI) | MetaboDesk |
|---|------|------------------------------|------------|
| 1 | **Model Yükleme** | `cobra.io.read_sbml_model()` | Open SBML butonu |
| 2 | **Model İstatistikleri** | `len(model.reactions)`, `len(model.metabolites)`, `len(model.genes)` | UI'daki sayaçlar |

### FAZA 2: Flux Analizleri
| # | Test | Geleneksel Araç | MetaboDesk |
|---|------|-----------------|------------|
| 3 | **FBA (Flux Balance Analysis)** | `model.optimize()` → `solution.objective_value`, `solution.fluxes` | FBA butonu |
| 4 | **pFBA (Parsimonious FBA)** | `cobra.flux_analysis.pfba(model)` | pFBA analiz tipi |
| 5 | **FVA (Flux Variability Analysis)** | `cobra.flux_analysis.flux_variability_analysis(model)` | FVA analiz tipi |
| 6 | **Flux Sampling** | `cobra.sampling.sample(model, n)` | Flux Sampling analiz tipi |

### FAZA 3: Gen/Reaksiyon Silme Analizleri
| # | Test | Geleneksel Araç | MetaboDesk |
|---|------|-----------------|------------|
| 7 | **Single Gene Deletion (SGD)** | `cobra.flux_analysis.single_gene_deletion(model)` | SGD analiz tipi |
| 8 | **Single Reaction Deletion (SRD)** | `cobra.flux_analysis.single_reaction_deletion(model)` | SRD analiz tipi |
| 9 | **Essential Reactions** | SRD sonuçlarından growth < %1 olanlar | Essential Reactions butonu |

### FAZA 4: Dual Değerler
| # | Test | Geleneksel Araç | MetaboDesk |
|---|------|-----------------|------------|
| 10 | **Shadow Prices** | `solution.shadow_prices` | Shadow Prices butonu |
| 11 | **Reduced Costs** | `solution.reduced_costs` | Reduced Costs butonu |

### FAZA 5: Kalite Testi
| # | Test | Geleneksel Araç | MetaboDesk |
|---|------|-----------------|------------|
| 12 | **Memote** | `memote run <model.xml>` (CLI) | Tools → Memote butonu |

---

## 🔄 Uygulama Sırası

### Adım 1: Test modeli indir
```python
import cobra
model = cobra.io.load_model("textbook")  # E. coli core
cobra.io.write_sbml_model(model, "e_coli_core.xml")
```

### Adım 2: Python scripti ile tüm analizleri çalıştır → CSV'lere kaydet
Bir Python scripti (`traditional_analysis.py`) yazılacak ve şunları yapacak:
1. Modeli yükle → istatistikleri kaydet
2. FBA çalıştır → objective value + tüm fluxları CSV'ye yaz
3. pFBA çalıştır → fluxları CSV'ye yaz
4. FVA çalıştır → min/max değerlerini CSV'ye yaz
5. SGD çalıştır → her gen silme sonucunu CSV'ye yaz
6. SRD çalıştır → her reaksiyon silme sonucunu CSV'ye yaz
7. Essential reactions listesini CSV'ye yaz
8. Shadow prices → CSV'ye yaz
9. Reduced costs → CSV'ye yaz
10. Flux sampling → istatistikleri CSV'ye yaz

### Adım 3: MetaboDesk ile aynı modeli aç ve aynı analizleri yap
- Her analiz sonucunu MetaboDesk'ten CSV olarak export et

### Adım 4: Karşılaştırma scripti yaz
Bir Python scripti (`compare_results.py`) yazılacak ve:
- Her iki tarafın CSV'lerini okuyacak
- Sayısal değerleri karşılaştıracak (tolerance: 1e-6)
- PASS/FAIL raporu üretecek
- Özet tablo oluşturacak

---

## 📊 Beklenen Çıktılar
```
validation/
├── VALIDATION_PLAN.md          ← Bu dosya
├── e_coli_core.xml             ← Test modeli
├── traditional_analysis.py     ← Geleneksel araç scripti
├── compare_results.py          ← Karşılaştırma scripti
├── traditional/                ← Geleneksel araç sonuçları
│   ├── model_info.json
│   ├── fba_fluxes.csv
│   ├── pfba_fluxes.csv
│   ├── fva_results.csv
│   ├── sgd_results.csv
│   ├── srd_results.csv
│   ├── essential_reactions.csv
│   ├── shadow_prices.csv
│   ├── reduced_costs.csv
│   └── sampling_stats.csv
├── metabodesk/                 ← MetaboDesk sonuçları (elle export edilecek)
│   └── (aynı dosya isimleri)
└── VALIDATION_REPORT.md        ← Karşılaştırma raporu
```

---

## ✅ Başarı Kriterleri
| Metrik | Kabul Edilebilir Fark |
|--------|----------------------|
| Objective Value (FBA) | ≤ 1e-6 |
| Flux değerleri | ≤ 1e-6 |
| FVA min/max | ≤ 1e-6 |
| SGD/SRD growth oranları | ≤ 1e-4 |
| Essential reaction listesi | Birebir aynı |
| Shadow prices | ≤ 1e-6 |
| Reduced costs | ≤ 1e-6 |

---

## 🗒️ Notlar
- Tüm analizler **aynı solver** (glpk) ile yapılacak
- Default medium/bounds kullanılacak (hiçbir değişiklik yapılmadan)
- Flux Sampling stokastik olduğu için ortalama ve standart sapma karşılaştırılacak
- Memote skoru manuel olarak karşılaştırılacak
