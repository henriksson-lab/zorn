mock5:
	cargo run --release -- --input /husky/henriksson/atrandi/v8_wgs_mock5_pta/countsketch_mat.csv --output umap_mock5.csv

s2:
	cargo run --release -- --input /husky/henriksson/atrandi/v8_wgs_saliva2_mda/countsketch_mat.csv --output umap_s2.csv

s2_sign:
	cargo run --release -- --input /husky/henriksson/atrandi/v8_wgs_saliva2_mda/countsketch_mat.csv --output umap_s2_sign.csv --transform sign

s2_log:
	cargo run --release -- --input /husky/henriksson/atrandi/v8_wgs_saliva2_mda/countsketch_mat.csv --output umap_s2_log.csv --transform log

s2_hamming:
	cargo run --release -- --input /husky/henriksson/atrandi/v8_wgs_saliva2_mda/countsketch_mat.csv --output umap_s2_hamming.csv --transform sign --metric hamming

s2_neighbors:
	cargo run --release -- --input /husky/henriksson/atrandi/v8_wgs_saliva2_mda/countsketch_mat.csv --output umap_s2_n30.csv --n-neighbors 30



s1:
	cargo run --release -- --input /husky/henriksson/atrandi/v8_wgs_saliva1_mda.shallower/countsketch_mat.csv --output umap_s1.csv

s1_sign:
	cargo run --release -- --input /husky/henriksson/atrandi/v8_wgs_saliva1_mda.shallower/countsketch_mat.csv --output umap_s1_sign.csv --transform sign

s1_log:
	cargo run --release -- --input /husky/henriksson/atrandi/v8_wgs_saliva1_mda.shallower/countsketch_mat.csv --output umap_s1_log.csv --transform log

s1_hamming:
	cargo run --release -- --input /husky/henriksson/atrandi/v8_wgs_saliva1_mda.shallower/countsketch_mat.csv --output umap_s1_hamming.csv --transform sign --metric hamming

s1_neighbors:
	cargo run --release -- --input /husky/henriksson/atrandi/v8_wgs_saliva1_mda.shallower/countsketch_mat.csv --output umap_s1_n30.csv --n-neighbors 30
