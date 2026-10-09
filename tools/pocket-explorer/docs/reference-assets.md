# Reference assets for a fresh checkout

The Crevice upgrade is committed as source, not a prebuilt application. The
current local preview already has its authoring inputs. A fresh checkout also
needs the generated organic anatomy, BodyParts3D atlas buffer and earlier
MakeHuman face buffer before `npm run build` or `npm run check`.

These reference buffers are 5,551,202 and 4,256,381 bytes. They exceed MARSE's
1 MiB repository limit and stay ignored; no size/privacy exception is added.
The smaller [reference-assets-lock.json](reference-assets-lock.json) records
their output hashes, public input identities, attribution and metadata without
bulk vertices. No local screenshots, native scenes or execution logs are added.

## BodyParts3D

Obtain the release 4.0 reduced archives from the official source linked in
[BODYPARTS3D_LICENSE.md](BODYPARTS3D_LICENSE.md). Put these public source files in
a local directory outside Git:

- `partof_BP3D_4.0_obj_99.zip`
- `isa_BP3D_4.0_obj_99.zip`
- `partof_element_parts.txt`

Check the archive SHA-256 identities against `records.atlas.metadata.archives`
in the lock, then run from the package directory:

```sh
python3 tools/build-atlas-context.py /path/to/public-atlas-directory
```

The builder makes no network requests. It writes ignored
`src/pocket-atlas-input.json`. Source element hashes are retained in both the
generated buffer and the lock. The atlas remains separate from authored pocket
geometry; its source units do not calibrate the crevice.

## Earlier graphical face dependency

The active Mouth/Face renderer uses the atlas. The older context module and its
tests still import the MakeHuman graphical face, so its buffer is required for
this source build. Use the pinned MakeHuman commit in `records.face.metadata`
and the asset licence recorded in [MAKEHUMAN_ASSET_LICENSE.md](MAKEHUMAN_ASSET_LICENSE.md).

Place the public inputs under the ignored `assets/face-reference/` directory as
`base.obj`, `african-mouth-open.target`, `asian-mouth-open.target`,
`caucasian-mouth-open.target` and `LICENSE.ASSETS.md`. Verify each input's SHA-256
against `records.face.metadata.inputs`. Run:

```sh
blender --background --factory-startup --disable-autoexec --python-exit-code 1 --python tools/build-face-reference.py
```

This writes ignored `src/pocket-face-input.json`. Registration and pose remain
authored display choices. The face is not a patient reconstruction.

## Finish the local build

Compare the SHA-256 of each generated JSON with its corresponding `sha256`
entry in the lock. A mismatch needs review; do not overwrite the lock to hide
it. Then run `npm run anatomy`, `npm run build` and `npm run check` as described
in the README. Blender 5.2.2 was used for the existing local outputs. The exact
fresh-checkout rebuild has not been replayed during commit preparation; the
current local build and 45 package tests have passed. Optional downloads are
operator actions; no source acquisition, upload or push occurs in this commit.
