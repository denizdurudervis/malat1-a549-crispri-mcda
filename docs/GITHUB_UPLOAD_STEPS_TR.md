# GitHub’a yükleme adımları

## Önerilen yöntem: GitHub Desktop

1. Bu ZIP paketini normal bir klasöre çıkar.
2. GitHub Desktop uygulamasını aç.
3. **File → Add local repository** seç.
4. Çıkardığın `MALAT1_A549_CRISPRi_MCDA_GitHub_v1` klasörünü göster.
5. Uygulama klasörün henüz Git deposu olmadığını söylerse **create a
   repository** seçeneğini kullan.
6. Önerilen depo adı:
   `malat1-a549-crispri-mcda`
7. İlk commit mesajı:
   `Initial reproducible research release`
8. **Publish repository** seç.
9. Makale yayımlanana kadar depoyu private tutmak mümkündür. Açık bilim
   paylaşımı planlandığında public yapılabilir.

## GitHub’a koyulmaması gerekenler

- `.venv-repro/`
- GRCh38 FASTA
- FlashFry bütün-genom indeks dosyaları
- büyük geçici dosyalar
- kişisel veya gizli veriler

Bunlar `.gitignore` ve provenance manifestleriyle yönetilir.

## Yayından sonra güncellenecek alanlar

- README içindeki repository URL
- `CITATION.cff`
- makale DOI’si
- hedef dergi ve sürüm bilgisi
- Zenodo arşiv DOI’si, oluşturulursa
