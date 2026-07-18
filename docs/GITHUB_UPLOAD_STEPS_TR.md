# GitHub’a yükleme adımları

## Önerilen yöntem: GitHub Desktop

1. Bu ZIP paketini normal bir klasöre çıkar.
2. GitHub Desktop uygulamasını aç.
3. **File → Add local repository** seç.
4. Çıkardığın `malat1-a549-crispri-mcda` klasörünü göster.
5. Uygulama klasörün henüz Git deposu olmadığını söylerse **create a
   repository** seçeneğini kullan.
6. Önerilen depo adı:
   `malat1-a549-crispri-mcda`
7. İlk commit mesajı:
   `Initial reproducible research release`
8. **Publish repository** seç.
9. Preprint yayımlanana kadar depoyu private tutmak mümkündür. Preprintteki
   yeniden üretilebilirlik bağlantısının çalışması için depo, preprint ile aynı
   gün public yapılmalı ve sürüm kalıcı bir arşive kaydedilmelidir.

## GitHub’a koyulmaması gerekenler

- `.venv-repro/`
- GRCh38 FASTA
- FlashFry bütün-genom indeks dosyaları
- büyük geçici dosyalar
- kişisel veya gizli veriler

Bunlar `.gitignore` ve provenance manifestleriyle yönetilir.

## Kalıcı tanımlayıcılar oluştuktan sonra güncellenecek alanlar

- makale DOI’si
- hedef dergi ve sürüm bilgisi
- Zenodo arşiv DOI’si
- `CITATION.cff` içindeki ilişkili makale ve arşiv alanları
