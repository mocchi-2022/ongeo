// Copyright (C) Mocchi (mocchi_2003@yahoo.co.jp)
// License: Boost Software License   See LICENSE.txt for the full license.
#include "ONGEO.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>

namespace {
const char *ParseString(const char *str, ON_String &parsed, const char delimiter, const char *default_v = 0){
	int idx_next, len;
	if (!str) goto ERROR_DEFAULT;
	if (str[0] == delimiter) goto ERROR_DEFAULT;
	if (std::sscanf(str, "%dH%n", &len, &idx_next) == 0) goto ERROR_DEFAULT;
	if (len) parsed.Append(str + idx_next, len);
	else if (default_v) parsed = default_v;
	const char *itr;
	if (delimiter != '\0') return (itr = std::strchr(str + idx_next + len, delimiter)) ? itr + 1 : 0;
	else return str + idx_next + len;
ERROR_DEFAULT:
	if (default_v) parsed = default_v;
	return 0;
}

const char *ParseInt(const char *str, int &parsed, const char delimiter, int default_v = 0){
	int idx_next;
	if (!str) goto ERROR_DEFAULT;
	if (str[0] == delimiter) goto ERROR_DEFAULT;
	if (std::sscanf(str, "%d%n", &parsed, &idx_next) == 0) goto ERROR_DEFAULT;
	const char *itr;
	if (delimiter != '\0') return (itr = std::strchr(str + idx_next, delimiter)) ? itr + 1 : 0;
	else return str + idx_next;
ERROR_DEFAULT:
	parsed = default_v;
	return 0;
}

const char *ParseDouble(const char *str, double &parsed, const char delimiter, double default_v = 0){
	int idx_next;
	if (!str) goto ERROR_DEFAULT;
	if (str[0] == delimiter) goto ERROR_DEFAULT;
	if (std::sscanf(str, "%lf%n", &parsed, &idx_next) == 0) goto ERROR_DEFAULT;
	const char *itr;
	if (delimiter != '\0') return (itr = std::strchr(str + idx_next, delimiter)) ? itr + 1 : 0;
	else return str + idx_next;
ERROR_DEFAULT:
	parsed = default_v;
	return 0;
}

void AppendIgesString(ON_String &dest, const ON_String &src, char delimit){
	if (src.Length() == 0){
		dest += delimit;
		return;
	}
	ON_String str;
	str.Format("%dH", src.Length()); str += src;
	dest += str;
	dest += delimit;
}

void AppendIgesInt(ON_String &dest, int src, char delimit){
	ON_String str;
	str.Format("%d", src);
	dest += str;
	dest += delimit;
}

void AppendIgesDouble(ON_String &dest, double src, char delimit){
	ON_String str;
	str.Format("%f", src);
	dest += str;
	dest += delimit;
}

bool ParseGlobalSection(const char *gs_str, ONGEO_IgesModel::GlobalSection &gs){
	const char *gs_iter = gs_str;
	if (*gs_iter == ',') gs.parm_delimiter = ",", ++gs_iter;
	else if (gs_iter[0] == '1' && gs_iter[1] == 'H' && gs_iter[2] == gs_iter[3]) gs.parm_delimiter = gs_iter[2], gs_iter += 4;
	else return false;
	if (*gs_iter == gs.parm_delimiter[0]) gs.record_delimiter = ";", ++gs_iter;
	else if (gs_iter[0] == '1' && gs_iter[1] == 'H' && gs_iter[3] == gs.parm_delimiter[0]) gs.record_delimiter = gs_iter[2], gs_iter += 4;
	else return false;

	gs_iter = ParseString(gs_iter, gs.product_id_sending,       gs.parm_delimiter[0]);
	gs_iter = ParseString(gs_iter, gs.file_name,                gs.parm_delimiter[0]);
	gs_iter = ParseString(gs_iter, gs.native_system_id,         gs.parm_delimiter[0]);
	gs_iter = ParseString(gs_iter, gs.preprocessor_version,     gs.parm_delimiter[0]);
	gs_iter = ParseInt(gs_iter, gs.num_of_binbits_for_intrep,   gs.parm_delimiter[0]);
	gs_iter = ParseInt(gs_iter, gs.max_pow_float,               gs.parm_delimiter[0]);
	gs_iter = ParseInt(gs_iter, gs.num_of_digits_float,         gs.parm_delimiter[0]);
	gs_iter = ParseInt(gs_iter, gs.max_pow_double,              gs.parm_delimiter[0]);
	gs_iter = ParseInt(gs_iter, gs.num_of_digits_double,        gs.parm_delimiter[0]);
	gs_iter = ParseString(gs_iter, gs.product_id_receiving,     gs.parm_delimiter[0], gs.product_id_sending.Array());
	gs_iter = ParseDouble(gs_iter, gs.model_scale,              gs.parm_delimiter[0], 1.0);
	gs_iter = ParseInt(gs_iter, gs.unit_flag,                   gs.parm_delimiter[0], 3);
	gs_iter = ParseString(gs_iter, gs.unit_name,                gs.parm_delimiter[0], "INCH");
	gs_iter = ParseInt(gs_iter, gs.max_num_lineweight_grad,     gs.parm_delimiter[0], 1);
	gs_iter = ParseDouble(gs_iter, gs.width_of_max_line_weight, gs.parm_delimiter[0]);
	gs_iter = ParseString(gs_iter, gs.timestamp_filegen,        gs.parm_delimiter[0]);
	gs_iter = ParseDouble(gs_iter, gs.min_resolution_in_unit,   gs.parm_delimiter[0]);
	gs_iter = ParseDouble(gs_iter, gs.approx_max_cod_in_unit,   gs.parm_delimiter[0], 0.0);
	gs_iter = ParseString(gs_iter, gs.name_author,              gs.parm_delimiter[0]);
	gs_iter = ParseString(gs_iter, gs.authors_org,              gs.parm_delimiter[0]);
	gs_iter = ParseInt(gs_iter, gs.flag_version,                gs.parm_delimiter[0], 3);
	gs_iter = ParseInt(gs_iter, gs.flag_draft_std,              gs.parm_delimiter[0], 0);
	gs_iter = ParseString(gs_iter, gs.timestamp_filemod,        gs.parm_delimiter[0]);
	gs_iter = ParseString(gs_iter, gs.desc_protocol,            gs.record_delimiter[0]);

	return true;
}

bool ParseDESection(const ON_String destr[2], ONGEO_IgesModel::DirectoryEntrySection &des){
	if (destr[0].Length() < 72 || destr[1].Length() < 72) return false;
	des.entity_type = std::atoi(destr[0].Left(8).Array());
	if (des.entity_type != std::atoi(destr[1].Left(8).Array())){
		des.entity_type = 0;
		return false;
	}

	des.param_data       = std::atoi(destr[0].Mid( 8, 8).Array());
	des.structure        = std::atoi(destr[0].Mid(16, 8).Array());
	des.line_font        = std::atoi(destr[0].Mid(24, 8).Array());
	des.level            = std::atoi(destr[0].Mid(32, 8).Array());
	des.view             = std::atoi(destr[0].Mid(40, 8).Array());
	des.trans_matrix     = std::atoi(destr[0].Mid(48, 8).Array());
	des.label_disp       = std::atoi(destr[0].Mid(56, 8).Array());
	des.stnum.blank_status  = std::atoi(destr[0].Mid(64, 2).Array());
	des.stnum.subord_ent_sw = std::atoi(destr[0].Mid(66, 2).Array());
	des.stnum.ent_use_flg   = std::atoi(destr[0].Mid(68, 2).Array());
	des.stnum.hierarchy     = std::atoi(destr[0].Mid(70, 2).Array());
	des.line_weight      = std::atoi(destr[1].Mid( 8, 8).Array());
	des.color_num        = std::atoi(destr[1].Mid(16, 8).Array());
	des.param_line_count = std::atoi(destr[1].Mid(24, 8).Array());
	des.form_num         = std::atoi(destr[1].Mid(32, 8).Array());
	des.reserved[0]      = std::atoi(destr[1].Mid(40, 8).Array());
	des.reserved[1]      = std::atoi(destr[1].Mid(48, 8).Array());
	::strncpy(des.ent_label, destr[1].Array() + 56, 8), des.ent_label[8] = '\0';
	des.ent_subscript    = std::atoi(destr[1].Mid(64, 8).Array());

	return true;
}
}

// 生成・破棄
ONGEO_IgesModel::ONGEO_IgesModel(){
}
ONGEO_IgesModel::ONGEO_IgesModel(const char *filename){
	Load(filename);
}
ONGEO_IgesModel::ONGEO_IgesModel(const ONGEO_IgesModel &rhs){
	ss = rhs.ss;
	gs = rhs.gs;
	des = rhs.des;
	ps = rhs.ps;
}
ONGEO_IgesModel &ONGEO_IgesModel::operator =(const ONGEO_IgesModel &rhs){
	ss = rhs.ss;
	gs = rhs.gs;
	des = rhs.des;
	ps = rhs.ps;
	return *this;
}
ONGEO_IgesModel::~ONGEO_IgesModel(){
}

bool ONGEO_IgesModel::Load(const char *filename){
	// 要件1 : 各行必ず80文字あること
	// 要件2 : 各行の73文字目がS,G,D,P,Tの順になっていること
	struct RAII{
		FILE *fp;
		RAII() : fp(0){}
		~RAII(){ if (fp) std::fclose(fp); }
	}raii;

	struct state_type{
		char mark;
		int count; // state毎の行数が入る。
		state_type(char mark_) : mark(mark_), count(0){}
	};
	state_type state[] = {state_type('S'), state_type('G'), state_type('D'), state_type('P'), state_type('T'), state_type('T')};
	state_type *current_state = &state[0];
	Clear();

	raii.fp = std::fopen(filename, "r");
	if (!raii.fp) return false;
	ON_SimpleArray<char> buf(82); buf.SetCount(buf.Capacity());
	int param_index = 0, param_rev_index = 0;

	ON_String global_section, de_section[2];

	for(;;){
		// 要件1チェック
		if (std::fread(buf.First(), 1, buf.Count(), raii.fp) < 80) return false;
		// 行末処理
		if (buf[81] != 0x0a) std::ungetc(buf[81], raii.fp);
		if (buf[80] != 0x0d && buf[80] != 0x0a) std::ungetc(buf[80], raii.fp);

		// 要件2チェック
		//  現在見ている行とcurrent_stateのmarkが異なる場合は、次のcurrent_stateに移す。
		//  移した先のcurrent_stateともmarkが異なる場合は要件2を満たしていないため、falseを返す。
		if (buf[72] != current_state->mark && buf[72] != (++current_state)->mark) return false;
		char cm = current_state->mark;
		if (cm == 'S'){
			ss.Append(&buf[0], 72);
		}else if (cm == 'G'){
			global_section.Append(&buf[0], 72);
		}else if (cm == 'D'){
			de_section[current_state->count & 1] = ON_String(&buf[0], 72);
			if (current_state->count & 1) ParseDESection(de_section, des.AppendNew());
		}else if (cm == 'P'){
			if (ps.Count() == 0){
				ps.SetCapacity(des.Count()), ps.SetCount(des.Count());
			}
			int param_rev_index_new = (std::sscanf(&buf[64],"%8uS", &param_rev_index_new), param_rev_index_new);
			ps[(param_rev_index_new-1)/2].Append(&buf[0], 64);
			param_rev_index = param_rev_index_new;
#if 0
			if (param_rev_index != param_rev_index_new){
//				pimpl->deptr.push_back(current_state->count+1);
				ps.AppendNew();
				param_rev_index = param_rev_index_new;
			}
			ps.Last()->Append(&buf[0], 64);
#endif
		}else if (cm == 'T'){
//			pimpl->deptr.push_back((current_state-1)->count+1);
			break;
		}
		++current_state->count;
	}

	ParseGlobalSection(global_section.Array(), gs);

	return (current_state->mark == 'T');
}

bool ONGEO_IgesModel::Save(const char *filename){
	// DictionaryEntryとParameterの数が異なる場合はエラーと認識する。
	if (des.Count() != ps.Count()) return false;
	// タイムスタンプの更新
	time_t tm = ::time(0);
	struct tm *now = ::localtime(&tm);
	gs.timestamp_filemod.Empty();
	gs.timestamp_filemod.Format("15H%04d%02d%02d.%02d%02d%02d",
		now->tm_year + 1900, now->tm_mon + 1, now->tm_mday,
		now->tm_hour, now->tm_min, now->tm_sec);

	// Parameter Sectionを整頓しなおす
	// 整頓により行数が変わるため、DirectoryEntryセクションのparam_dataも修正する。
	if (des.Count()) des[0].param_data = 1;
	for (int j = 0; j < ps.Count(); ++j){
		ON_String str = ps[j];
		str.Replace(" ", "");
		const char *p = str.Array();
		const char *pe = p + str.Length();

		ON_String str_mod;
		const char *pi = p, *pl = p;
		while(pi < pe){
			const char *pq = ::strchr(pi, gs.parm_delimiter[0]);
			pq = pq ? pq + 1 : pe;
			size_t llen = pq - pl, dif = pq - pi;
			while (dif > 64){
				str_mod.Append(pi, 64);
				pi += 64;
				dif -= 64;
				llen -= 64;
			}
			if (llen > 64){
				llen = pi - pl;
				str_mod.Append(pl, static_cast<int>(llen));
				if (llen < 64) str_mod += ON_String(' ', 64 - llen);
				pl = pi;
			}
			pi += dif;
		}
		str_mod += pl;
		ps[j] = str_mod;
		des[j].param_line_count = (str_mod.Length() + 63) / 64;
		if (j > 0) des[j].param_data = des[j-1].param_data + des[j-1].param_line_count;
		str_mod.Empty();
	}

	FILE *fp = std::fopen(filename, "w");
	if (!fp) return false;

	// Globalセクションの文字列化
	ON_String gss;
	if (gs.parm_delimiter == ",") gss += ",";
	else AppendIgesString(gss, gs.parm_delimiter, gs.parm_delimiter[0]);
	if (gs.record_delimiter == ";") gss += gs.parm_delimiter;
	else AppendIgesString(gss, gs.record_delimiter, gs.parm_delimiter[0]);
	AppendIgesString(gss, gs.product_id_sending, gs.parm_delimiter[0]);
	AppendIgesString(gss, gs.file_name, gs.parm_delimiter[0]);
	AppendIgesString(gss, gs.native_system_id, gs.parm_delimiter[0]);
	AppendIgesString(gss, gs.preprocessor_version, gs.parm_delimiter[0]);
	AppendIgesInt   (gss, gs.num_of_binbits_for_intrep, gs.parm_delimiter[0]);
	AppendIgesInt   (gss, gs.max_pow_float, gs.parm_delimiter[0]);
	AppendIgesInt   (gss, gs.num_of_digits_float, gs.parm_delimiter[0]);
	AppendIgesInt   (gss, gs.max_pow_double, gs.parm_delimiter[0]);
	AppendIgesInt   (gss, gs.num_of_digits_double, gs.parm_delimiter[0]);
	AppendIgesString(gss, gs.product_id_receiving, gs.parm_delimiter[0]);
	AppendIgesDouble(gss, gs.model_scale, gs.parm_delimiter[0]);
	AppendIgesInt   (gss, gs.unit_flag, gs.parm_delimiter[0]);
	AppendIgesString(gss, gs.unit_name, gs.parm_delimiter[0]);
	AppendIgesInt   (gss, gs.max_num_lineweight_grad, gs.parm_delimiter[0]);
	AppendIgesDouble(gss, gs.width_of_max_line_weight, gs.parm_delimiter[0]);
	AppendIgesString(gss, gs.timestamp_filegen, gs.parm_delimiter[0]);
	AppendIgesDouble(gss, gs.min_resolution_in_unit, gs.parm_delimiter[0]);
	AppendIgesDouble(gss, gs.approx_max_cod_in_unit, gs.parm_delimiter[0]);
	AppendIgesString(gss, gs.name_author, gs.parm_delimiter[0]);
	AppendIgesString(gss, gs.authors_org, gs.parm_delimiter[0]);
	AppendIgesInt   (gss, gs.flag_version, gs.parm_delimiter[0]);
	AppendIgesInt   (gss, gs.flag_draft_std, gs.parm_delimiter[0]);
	AppendIgesString(gss, gs.timestamp_filemod, gs.parm_delimiter[0]);
	AppendIgesString(gss, gs.desc_protocol, gs.record_delimiter[0]);

	// Startセクション、Globalセクション
	int lnums[] = {0, 0, 0, 0};
	char mark[] = {'S', 'G'};
	for (int i = 0; i < 2; ++i){
		ON_String &s = (i == 0) ? ss : gss;
		int &lnum = lnums[i], slen = s.Length();
		while(slen){
			++lnum;
			if (slen >= 72){
				fprintf(fp, "%s%c%7d\n", s.Mid((lnum - 1) * 72, 72), mark[i], lnum);
				slen -= 72;
			}else{
				fprintf(fp, "%s%s%c%7d\n", s.Mid((lnum - 1) * 72), ON_String(' ', 72 - slen), mark[i], lnum);
				slen = 0;
			}
		}
	}
	// DirectoryEntryセクション
	for (int i = 0; i < des.Count(); ++i){
		DirectoryEntrySection &d = des[i];
		std::fprintf(fp, "%8d%8d%8d%8d%8d%8d%8d%8d%02d%02d%02d%02dD%7d\n",
			d.entity_type, d.param_data, d.structure, d.line_font, d.level, d.view, d.trans_matrix, d.label_disp,
			d.stnum.blank_status, d.stnum.subord_ent_sw, d.stnum.ent_use_flg, d.stnum.hierarchy, ++lnums[2]);
		std::fprintf(fp, "%8d%8d%8d%8d%8d%8d%8d%8s%8dD%7d\n",
			d.entity_type, d.line_weight, d.color_num, d.param_line_count, d.form_num, d.reserved[0], d.reserved[1],
			d.ent_label, d.ent_subscript, ++lnums[2]);
	}
	// Parameterセクション
	for (int i = 0; i < ps.Count(); ++i){
		ON_String &s = ps[i];
		int slen = s.Length(), &lnum = lnums[3], ldes = i * 2 + 1;
		int sindex = 0;
		while(slen){
			++lnum;
			if (slen >= 64){
				std::fprintf(fp, "%s %7dP%7d\n", s.Mid(sindex, 64), ldes, lnum);
				slen -= 64;
				sindex += 64;
			}else{
				std::fprintf(fp, "%s%s %7dP%7d\n", s.Mid(sindex), ON_String(' ', 64 - slen), ldes, lnum);
				slen = 0;
			}
		}
	}
	// Terminalセクション
	std::fprintf(fp, "S%7dG%7dD%7dP%7d%sT      1\n", lnums[0], lnums[1], lnums[2], lnums[3], ON_String(' ', 40));

	::fclose(fp);
	return true;
}
void ONGEO_IgesModel::Clear(){
	ss.Empty();
	gs = GlobalSection();
	des.Empty();
	ps.Empty();
}

bool ONGEO_IgesModel::SetEntityClearOut(int index){
	if (index < 0 || index >= des.Count()) return false;
	des[index] = DirectoryEntrySection();
	ps[index].Empty();
	return true;
}



ONGEO_IgesModel *ONGEO_NewIgesModel(){
	return new ONGEO_IgesModel();
}
ONGEO_IgesModel *ONGEO_NewIgesModel(const char *filename){
	return new ONGEO_IgesModel(filename);
}
void ONGEO_DeleteIgesModel(ONGEO_IgesModel *oim){
	delete oim;
}
bool ONGEO_IgesModel_Load(ONGEO_IgesModel *oim, const char *filename){
	if (!oim) return false;
	return oim->Load(filename);
}
bool ONGEO_IgesModel_Save(ONGEO_IgesModel *oim, const char *filename){
	if (!oim) return false;
	return oim->Save(filename);
}

void ONGEO_IgesModel_Clear(ONGEO_IgesModel *oim){
	if (oim) oim->Clear();
}

void ONGEO_IgesModel_SetEntityClearOut(ONGEO_IgesModel *oim, int index){
	if (oim) oim->SetEntityClearOut(index);
}
