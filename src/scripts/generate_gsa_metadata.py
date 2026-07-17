import pandas as pd
import os
import hashlib
import argparse


def parse_args():
    """
    解析命令行参数

    返回:
        包含fq.gz文件存储路径等参数的命名空间
    """
    parser = argparse.ArgumentParser(
        description='生成 GSA-Human HRA 元数据提交表格'
    )
    parser.add_argument(
        '--fq-dir',
        type=str,
        default='data/fastq',
        help='fq.gz 文件存储目录路径（默认: data/fastq）'
    )
    parser.add_argument(
        '--plasma-exp',
        type=str,
        default='data/plasma_exp.xlsx',
        help='plasma_exp.xlsx 文件路径（默认: data/plasma_exp.xlsx）'
    )
    parser.add_argument(
        '--samples-v4',
        type=str,
        default='data/ovarian_fq_series_samples_v4.xlsx',
        help='ovarian_fq_series_samples_v4.xlsx 文件路径（默认: data/ovarian_fq_series_samples_v4.xlsx）'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='GSA_HRA_Metadata_Submission.xlsx',
        help='输出文件路径（默认: GSA_HRA_Metadata_Submission.xlsx）'
    )
    return parser.parse_args()


def calculate_md5(file_path):
    """
    计算文件的MD5值

    参数:
        file_path: 文件路径

    返回:
        MD5字符串，如果文件不存在返回空字符串
    """
    if not os.path.exists(file_path):
        return ''

    hash_md5 = hashlib.md5()
    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def map_to_doid(pathology):
    """
    将 pathology 自由文本映射到 DO (Disease Ontology) 术语

    参数:
        pathology: 原始 pathology 文本

    返回:
        DOID 术语，如 'DOID:0050'
    """
    if not pathology:
        return ''

    pathology_lower = pathology.lower()

    # DO (Disease Ontology) 卵巢癌相关术语映射
    # 参考: https://disease-ontology.org/
    # 顺序重要：更具体的术语放在前面
    do_mapping = [
        ('borderline', 'DOID:162'),  # ovarian borderline tumor (最优先)
        ('clear cell', 'DOID:5998'),  # ovarian clear cell adenocarcinoma
        ('endometrioid', 'DOID:5615'),  # ovarian endometrioid adenocarcinoma
        ('mucinous', 'DOID:4960'),  # ovarian mucinous adenocarcinoma
        ('low-grade', 'DOID:4961'),  # low-grade serous adenocarcinoma
        ('serous', 'DOID:4961'),  # ovarian serous adenocarcinoma
        ('malignant', 'DOID:2394'),  # ovarian carcinoma
        ('advanced', 'DOID:0050'),  # ovarian cancer
        ('high-grade', 'DOID:2394'),  # high-grade serous carcinoma
    ]

    for key, doid in do_mapping:
        if key in pathology_lower:
            return doid

    # 默认返回卵巢癌通用术语
    return 'DOID:0050'  # ovarian cancer


def generate_individual(df_plasma, df_plasma_v4):
    """
    生成 Individual sheet 数据

    使用 GSA-Human HRA 模板的完整列名
    ID 使用 case_id（加 'D' 前缀）
    """
    # 模板完整列名
    columns = [
        '*ID', '*individual_identifier', '*gender', 'existed_individual_accession',
        '*experiment_group', 'disease_name', 'TCM_syndrome_name', 'height (m)',
        'body_weight (kg)', 'body_mass_index(BMI)', 'province', 'occupation',
        'birthplace', 'race', 'vital_status', 'age', 'cause_of_death',
        'days_to_birth', 'days_to_death', 'year_of_birth', 'year_of_death',
        'primary_diagnosis', 'age_at_diagnosis', 'diagnosis_is_primary_disease',
        'morphology', 'tissue_or_organ_of_origin', 'ajcc_clinical_stage',
        'ajcc_clinical_treatment', 'classification_of_tumor', 'metastasis_at_diagnosis',
        'prior_treatment', 'residual_disease', 'tumor_grade', 'is_drinker', 'is_smoker',
        'relationship_type', 'relationship_age_at_diagnosis', 'relationship_primary_diagnosis',
        'relative_deceased', 'days_to_follow_up', 'days_to_progression',
        'days_to_recurrence', 'history_of_tumor', 'progression_or_recurrence',
        'progression_or_recurrence_type', 'chemo_concurrent_to_radiation',
        'drug_category', 'initial_disease_status', 'therapeutic_agents',
        'treatment_dose', 'treatment_dose_units', 'treatment_effect',
        'treatment_frequency', 'treatment_outcome', 'tumor_burden',
        'tumor_largest_dimension_diameter', 'tumor_length_measurement',
        'tumor_shape', 'tumor_thickness', 'custom_attributes_1',
    ]

    individual_data = []

    for _, row in df_plasma_v4.iterrows():
        case_id = str(row['case_id'])

        pathology = row['pathology'] if 'pathology' in row and pd.notna(row['pathology']) else ''
        subtype = row['subtype'] if 'subtype' in row and pd.notna(row['subtype']) else ''
        category = row['category'] if 'category' in row and pd.notna(row['category']) else ''

        data = {
            '*ID': 'D' + case_id,
            '*individual_identifier': case_id,
            '*gender': 'female',
            'existed_individual_accession': '',
            '*experiment_group': 'Case',
            'disease_name': map_to_doid(pathology),  # DOID 格式
            'TCM_syndrome_name': '',
            'height (m)': '',
            'body_weight (kg)': '',
            'body_mass_index(BMI)': '',
            'province': '',
            'occupation': '',
            'birthplace': '',
            'race': '',
            'vital_status': '',
            'age': '',
            'cause_of_death': '',
            'days_to_birth': '',
            'days_to_death': '',
            'year_of_birth': '',
            'year_of_death': '',
            'primary_diagnosis': pathology,  # 原始病理文本
            'age_at_diagnosis': '',
            'diagnosis_is_primary_disease': '',
            'morphology': '',
            'tissue_or_organ_of_origin': '',
            'ajcc_clinical_stage': '',
            'ajcc_clinical_treatment': '',
            'classification_of_tumor': category,
            'metastasis_at_diagnosis': '',
            'prior_treatment': '',
            'residual_disease': '',
            'tumor_grade': 'G3' if category == 'malignant' else '',
            'is_drinker': '',
            'is_smoker': '',
            'relationship_type': '',
            'relationship_age_at_diagnosis': '',
            'relationship_primary_diagnosis': '',
            'relative_deceased': '',
            'days_to_follow_up': '',
            'days_to_progression': '',
            'days_to_recurrence': '',
            'history_of_tumor': '',
            'progression_or_recurrence': '',
            'progression_or_recurrence_type': '',
            'chemo_concurrent_to_radiation': '',
            'drug_category': '',
            'initial_disease_status': '',
            'therapeutic_agents': '',
            'treatment_dose': '',
            'treatment_dose_units': '',
            'treatment_effect': '',
            'treatment_frequency': '',
            'treatment_outcome': '',
            'tumor_burden': '',
            'tumor_largest_dimension_diameter': '',
            'tumor_length_measurement': '',
            'tumor_shape': '',
            'tumor_thickness': '',
            'custom_attributes_1': '',
        }

        individual_data.append(data)

    df = pd.DataFrame(individual_data)
    for col in columns:
        if col not in df.columns:
            df[col] = ''

    return df[columns]

def generate_sample(df_plasma, df_plasma_v4):
    """
    生成 Sample sheet 数据

    使用 GSA-Human HRA 模板的完整列名
    ID 使用 case_id（加 'S' 前缀）
    """
    # 模板完整列名
    columns = [
        '*ID', '*sample_name', '*tissue', '*individual_ID', 'age', 'age_unit',
        '*public_description', '*sample_title', '*collection_date', '*biomaterial_provider',
        '*geographic_location', 'culture_collection', 'karyotype', 'phenotype', 'population',
        'custom_attributes_1',
    ]

    sample_data = []

    for _, row in df_plasma.iterrows():
        case_id = str(row['case_id'])

        sample_name = row['id'] if 'id' in row and pd.notna(row['id']) else case_id

        collection_date = ''
        if 'datetime' in row and pd.notna(row['datetime']):
            dt_str = str(row['datetime'])
            if len(dt_str) >= 8:
                collection_date = dt_str[:4] + '-' + dt_str[4:6] + '-' + dt_str[6:8]

        data = {
            '*ID': 'S' + case_id,
            '*sample_name': sample_name,
            '*tissue': 'Blood',
            '*individual_ID': 'D' + case_id,
            'age': '',
            'age_unit': '',
            '*public_description': 'Cell-free DNA from plasma of ovarian cancer patient',
            '*sample_title': 'Plasma sample ' + case_id,
            '*collection_date': collection_date,
            '*biomaterial_provider': 'Beijing hospital',
            '*geographic_location': 'Beijing, China',
            'culture_collection': '',
            'karyotype': '',
            'phenotype': '',
            'population': 'Chinese',
            'custom_attributes_1': '',
        }

        sample_data.append(data)

    df = pd.DataFrame(sample_data)
    for col in columns:
        if col not in df.columns:
            df[col] = ''

    return df[columns]

def generate_experiment(sample_df):
    """
    生成 Experiment sheet 数据

    使用 GSA-Human HRA 模板的完整列名
    ID 使用 sample ID（加 'E' 前缀）
    """
    # 模板完整列名
    columns = [
        '*ID', '*Experiment title', '*sample_ID', '*Platform',
        '*Library Construction / Experimental Design', 'Library name',
        '*Strategy', '*Source', '*Selection', '*Layout',
        '*Read length for mate1(bp)', 'Read length for mate 2(bp)',
        'Insert size (bp)', 'Nominal size (bp)', 'Nominal standard deviation (bp)',
        'Planned number of cycles',
    ]

    experiment_data = []

    for _, row in sample_df.iterrows():
        sample_id = row['*ID']
        sample_name = row['*sample_name']

        data = {
            '*ID': 'E' + sample_id[1:] if sample_id.startswith('S') else sample_id,
            '*Experiment title': 'Methylation sequencing of ' + sample_name,
            '*sample_ID': sample_id,
            '*Platform': 'Illumina NovaSeq 6000',
            '*Library Construction / Experimental Design': 'Bisulfite sequencing',
            'Library name': '',
            '*Strategy': 'WGS',
            '*Source': 'GENOMIC',
            '*Selection': 'RANDOM',
            '*Layout': 'PAIRED',
            '*Read length for mate1(bp)': 150,
            'Read length for mate 2(bp)': 150,
            'Insert size (bp)': '',
            'Nominal size (bp)': '',
            'Nominal standard deviation (bp)': '',
            'Planned number of cycles': '',
        }

        experiment_data.append(data)

    df = pd.DataFrame(experiment_data)
    for col in columns:
        if col not in df.columns:
            df[col] = ''

    return df[columns]

def generate_run(experiment_df, df_plasma_v4, data_dir='data/fastq'):
    """
    生成 Run sheet 数据

    使用 GSA-Human HRA 模板的完整列名
    ID 使用 experiment ID（加 'R' 前缀）
    文件名使用 plasma_v4 中的 fq_R1/fq_R2
    计算文件MD5值
    """
    # 模板完整列名
    columns = [
        '*ID', '*Run title', '*Experiment ID', '*Run data file type',
        '*File name 1', '*MD5 checksum 1', 'File name 2', 'MD5 checksum 2',
        'Reference file name', 'MD5 for reference file',
        'Assembly Name or Accession', 'Assembly Accession URL',
    ]

    run_data = []

    # case_id 映射（用于查找文件名）
    # Sample ID 格式是 S640849，需要提取 640849
    sample_id_to_case_id = {}
    for _, row in experiment_df.iterrows():
        sample_id = row['*sample_ID']
        if sample_id.startswith('S'):
            case_id = sample_id[1:]
            sample_id_to_case_id[case_id] = case_id

    case_id_to_fq = {}
    for _, row in df_plasma_v4.iterrows():
        case_id = str(row['case_id'])
        case_id_to_fq[case_id] = {
            'fq_R1': row['fq_R1'] if 'fq_R1' in row and pd.notna(row['fq_R1']) else '',
            'fq_R2': row['fq_R2'] if 'fq_R2' in row and pd.notna(row['fq_R2']) else '',
        }

    for _, row in experiment_df.iterrows():
        exp_id = row['*ID']
        sample_id = row['*sample_ID']

        # 从 sample_id 提取 case_id
        case_id = sample_id[1:] if sample_id.startswith('S') else sample_id

        fq_info = case_id_to_fq.get(case_id, {})

        fq_r1 = fq_info.get('fq_R1', '')
        fq_r2 = fq_info.get('fq_R2', '')

        md5_1 = ''
        md5_2 = ''

        if fq_r1:
            fq_path_1 = os.path.join(data_dir, fq_r1)
            md5_1 = calculate_md5(fq_path_1)
            if md5_1:
                print('[OK] {} -> {}'.format(fq_r1, md5_1))
            else:
                print('[MISSING] {}'.format(fq_path_1))

        if fq_r2:
            fq_path_2 = os.path.join(data_dir, fq_r2)
            md5_2 = calculate_md5(fq_path_2)
            if md5_2:
                print('[OK] {} -> {}'.format(fq_r2, md5_2))
            else:
                print('[MISSING] {}'.format(fq_path_2))

        data = {
            '*ID': 'R' + exp_id[1:] if exp_id.startswith('E') else exp_id,
            '*Run title': 'Read data for ' + sample_id,
            '*Experiment ID': exp_id,
            '*Run data file type': 'fastq',
            '*File name 1': fq_r1,
            '*MD5 checksum 1': md5_1,
            'File name 2': fq_r2,
            'MD5 checksum 2': md5_2,
            'Reference file name': '',
            'MD5 for reference file': '',
            'Assembly Name or Accession': 'GRCh38',
            'Assembly Accession URL': 'https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/',
        }

        run_data.append(data)

    df = pd.DataFrame(run_data)
    for col in columns:
        if col not in df.columns:
            df[col] = ''

    return df[columns]

def main():
    """
    主函数：生成 GSA-Human HRA 元数据表格（plasma 样本交集）

    支持通过命令行参数指定 fq.gz 文件存储路径、输入文件和输出文件路径
    """
    args = parse_args()

    FQ_DIR = args.fq_dir

    if not os.path.exists(FQ_DIR):
        print('Warning: Data directory does not exist: {}'.format(FQ_DIR))
        print('MD5 calculation will be skipped.')

    df_plasma = pd.read_excel(args.plasma_exp, sheet_name='plasma')
    df_plasma_v4 = pd.read_excel(args.samples_v4, sheet_name='plasma')

    print('plasma_exp: {} samples'.format(len(df_plasma)))
    print('plasma_v4: {} samples'.format(len(df_plasma_v4)))

    plasma_case_ids = set(df_plasma['case_id'].astype(str).unique())
    plasma_v4_case_ids = set(df_plasma_v4['case_id'].astype(str).unique())
    common_case_ids = plasma_case_ids.intersection(plasma_v4_case_ids)
    print('Common case_ids (intersection): {}'.format(len(common_case_ids)))

    df_plasma = df_plasma[df_plasma['case_id'].astype(str).isin(common_case_ids)]
    df_plasma_v4 = df_plasma_v4[df_plasma_v4['case_id'].astype(str).isin(common_case_ids)]

    print('Filtered plasma_exp: {} samples'.format(len(df_plasma)))
    print('Filtered plasma_v4: {} samples'.format(len(df_plasma_v4)))

    individual_df = generate_individual(df_plasma, df_plasma_v4)
    print('\nGenerated Individual: {} individuals'.format(len(individual_df)))

    sample_df = generate_sample(df_plasma, df_plasma_v4)
    print('Generated Sample: {} samples'.format(len(sample_df)))

    experiment_df = generate_experiment(sample_df)
    print('Generated Experiment: {} experiments'.format(len(experiment_df)))

    run_df = generate_run(experiment_df, df_plasma_v4, FQ_DIR)
    print('Generated Run: {} runs'.format(len(run_df)))

    md5_count = (run_df['*MD5 checksum 1'] != '').sum() + (run_df['MD5 checksum 2'] != '').sum()
    print('Calculated MD5 for {} files'.format(md5_count))

    output_file = args.output
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        individual_df.to_excel(writer, sheet_name='Individual', index=False)
        sample_df.to_excel(writer, sheet_name='Sample', index=False)
        experiment_df.to_excel(writer, sheet_name='Experiment', index=False)
        run_df.to_excel(writer, sheet_name='Run', index=False)

    print('\nMetadata file saved to: {}'.format(output_file))


if __name__ == '__main__':
    main()
