# BodyParts3D oral and craniofacial context

BodyParts3D, © The Database Center for Life Science licensed under CC Attribution 4.0 International.

Official source: https://dbarchive.biosciencedbc.jp/en/bodyparts3d/download.html

Official license, updated 2025-02-27: https://dbarchive.biosciencedbc.jp/en/bodyparts3d/lic.html

License text: https://creativecommons.org/licenses/by/4.0/legalcode

Data: release 4.0 (2013), official IS-A and PART-OF OBJ archives with the
archive's stated 99% polygon reduction. The old OBJ headers still name CC BY-SA
2.1 Japan. Those original bytes and hashes are preserved locally; the official
archive now identifies this database's license as CC BY 4.0.

The derived viewer subset selects 63 source elements, triangulates polygon
faces and packs original positions as float32 and indices as uint32. Source
FMA concepts, element identifiers, archive identity and original OBJ SHA-256
are retained in `src/pocket-atlas-input.json`. No source code was imported.

Original atlas coordinates remain common to every source part. Display uses a
rigid basis change: (x, z−1460, −y−140). The source OBJ headers state millimeters.
An authored lower-jaw display opening rotates the mandible, lower teeth,
lower gingiva and tongue together by 0.45 radians about (0,40,−45) in viewer
coordinates. This is not measured TMJ movement. Other soft-tissue deformation
has not been modeled. The atlas is an adult male reference, not an average
adult mouth or patient reconstruction. Expert anatomical review is pending.

The selected atlas FDI 36 is identified from the source name “left lower first
secondary molar tooth,” FMA55704, element FJ1254. The detailed periodontal
pocket is an independent authored specimen. The two have no measured geometric
registration. Whole atlas teeth are single surfaces, not separate enamel,
dentin and pulp segmentations. No species locations or disease presets are
assigned to atlas tissues.

MakeHuman facial assets from the earlier graphical experiment are preserved
locally with their CC0 receipt but are not used by the active Face/Mouth views.
