//! Dense-region filter. OWNER: P1.
//!
//! Mirrors src/combine_insertions_region_filter.py `filter_dense_regions` (called with
//! bin_range=100, ins_cutoff=4 from combine_insertions.py:139-141). SPEC.md §3.4.

use crate::model::Insertion;
use rustc_hash::FxHashMap;
use rustc_hash::FxHashSet;

/// Returns (kept in input order, n_removed, n_hot_bins).
///
/// pos = right_pos if present else left_pos. bin = `int(pos / bin_range) * bin_range` (python
/// true division then truncation; positions are >= 0). Hot bins: (contig, bin) with count >=
/// cutoff. An insertion is removed iff some hot bin b on its contig satisfies
/// `pos > b - bin_range/2 and pos < b + bin_range*1.5` (float comparisons; with bin_range=100
/// these are exact integers 50 / 150). Candidate bins: b from
/// `((pos - int(bin_range*1.5)) // bin_range) * bin_range` up to (exclusive)
/// `pos + int(bin_range/2) + bin_range`, step bin_range.
pub fn filter_dense_regions(insertions: Vec<Insertion>, bin_range: i64, ins_cutoff: usize) -> (Vec<Insertion>, usize, usize) {
    fn pos_of(i: &Insertion) -> i64 {
        i.right_pos.or(i.left_pos).expect("filter_dense_regions: insertion without a coordinate")
    }
    let br = bin_range as f64;
    let mut counts: FxHashMap<(u32, i64), usize> = FxHashMap::default();
    for i in &insertions {
        let b = ((pos_of(i) as f64 / br) as i64) * bin_range;
        *counts.entry((i.contig, b)).or_insert(0) += 1;
    }
    let hot: FxHashSet<(u32, i64)> = counts.into_iter().filter(|&(_, n)| n >= ins_cutoff).map(|(k, _)| k).collect();
    let n_hot = hot.len();
    let lead = (br * 1.5) as i64; // int(bin_range * 1.5)
    let half = (br / 2.0) as i64; // int(bin_range / 2)
    let mut kept = Vec::with_capacity(insertions.len());
    let mut removed = 0usize;
    for i in insertions {
        let pos = pos_of(&i);
        let base = (pos - lead).div_euclid(bin_range) * bin_range;
        let stop = pos + half + bin_range;
        let mut hit = false;
        let mut bin = base;
        while bin < stop {
            if hot.contains(&(i.contig, bin)) && (pos as f64) > bin as f64 - br / 2.0 && (pos as f64) < bin as f64 + br * 1.5 {
                hit = true;
                break;
            }
            bin += bin_range;
        }
        if hit {
            removed += 1;
        } else {
            kept.push(i);
        }
    }
    (kept, removed, n_hot)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{InsType, Interner, Tok};

    fn mk(contigs: &Interner, c: &str, left: i64, right: Option<i64>) -> Insertion {
        Insertion {
            uid: 0,
            contig: contigs.intern(c),
            name_start: Tok::pos(left),
            name_end: Tok::pos(right.unwrap_or(0)),
            ty: InsType::FullInfo,
            open_side: None,
            left_clipped: None,
            left_aligned: None,
            left_pos: Some(left),
            right_clipped: None,
            right_aligned: None,
            right_pos: right,
            files: vec![],
            member_loci: vec![],
            member_sides: None,
        }
    }

    /// expected values from python `filter_dense_regions(ins, 100, cut)`
    fn check_rf(ins: &[(&str, i64, Option<i64>)], cut: usize, kept: &[usize], removed: usize, hot: usize) {
        let contigs = Interner::new();
        let v: Vec<Insertion> = ins.iter().enumerate().map(|(n, (c, l, r))| {
            let mut i = mk(&contigs, c, *l, *r);
            i.uid = n as u32;
            i
        }).collect();
        let (k, rem, nh) = filter_dense_regions(v, 100, cut);
        assert_eq!(k.iter().map(|i| i.uid as usize).collect::<Vec<_>>(), kept);
        assert_eq!((rem, nh), (removed, hot));
    }

    #[test]
    fn dense_regions_match_python() {
        check_rf(&[("chr1", 1247, None), ("chr2", 1247, None), ("chr2", 1593, Some(1603)), ("chr2", 1487, None), ("chr1", 1686, None), ("chr1", 1463, None), ("chr1", 1038, None), ("chr1", 1192, Some(1202)), ("chr1", 1339, None), ("chr2", 1654, Some(1664)), ("chr2", 1075, None), ("chr2", 1250, None), ("chr1", 1302, Some(1312)), ("chr1", 1515, Some(1525)), ("chr1", 1591, Some(1601)), ("chr1", 1663, None), ("chr1", 1291, None), ("chr2", 1038, Some(1048)), ("chr1", 1512, Some(1522)), ("chr1", 1675, Some(1685)), ("chr2", 1667, Some(1677)), ("chr1", 1468, Some(1478)), ("chr1", 1667, Some(1677)), ("chr1", 1407, Some(1417)), ("chr2", 1363, Some(1373)), ("chr2", 1677, None), ("chr1", 1028, Some(1038)), ("chr2", 1060, Some(1070)), ("chr1", 1058, None), ("chr2", 1633, Some(1643)), ("chr2", 1249, Some(1259)), ("chr2", 1576, None), ("chr2", 1064, None), ("chr2", 1119, Some(1129)), ("chr2", 1437, Some(1447)), ("chr1", 1189, None), ("chr1", 1066, Some(1076)), ("chr2", 1235, Some(1245)), ("chr2", 1346, Some(1356)), ("chr2", 1416, None), ("chr2", 1165, None), ("chr2", 1173, Some(1183)), ("chr1", 1556, None), ("chr1", 1652, None), ("chr2", 1109, Some(1119)), ("chr2", 1094, None), ("chr1", 1146, Some(1156)), ("chr2", 1086, None), ("chr2", 1572, None), ("chr2", 1659, None), ("chr1", 1318, None), ("chr2", 1335, Some(1345)), ("chr2", 1522, None), ("chr1", 1624, Some(1634)), ("chr1", 1505, Some(1515)), ("chr2", 1680, None), ("chr1", 1534, Some(1544)), ("chr1", 1643, None), ("chr2", 1519, Some(1529)), ("chr2", 1547, Some(1557)), ("chr2", 1032, None), ("chr2", 1328, None), ("chr1", 1462, None), ("chr2", 1455, Some(1465)), ("chr2", 1659, Some(1669)), ("chr1", 1442, None), ("chr2", 1640, Some(1650)), ("chr1", 1279, None), ("chr2", 1359, Some(1369)), ("chr2", 1622, None), ("chr2", 1224, Some(1234)), ("chr2", 1111, None), ("chr1", 1489, Some(1499)), ("chr1", 1276, None), ("chr1", 1134, None), ("chr2", 1439, Some(1449)), ("chr1", 1311, Some(1321)), ("chr1", 1151, Some(1161)), ("chr1", 1690, None), ("chr2", 1037, None)], 2, &[], 80, 14);
        check_rf(&[("chr2", 1424, Some(1434)), ("chr2", 1511, Some(1521)), ("chr2", 1450, None), ("chr2", 1010, Some(1020)), ("chr1", 1388, None), ("chr2", 1116, None), ("chr2", 1129, None), ("chr1", 1201, Some(1211)), ("chr1", 1290, Some(1300)), ("chr1", 1228, Some(1238)), ("chr2", 1114, Some(1124)), ("chr2", 1528, None), ("chr1", 1512, Some(1522)), ("chr1", 1207, Some(1217)), ("chr2", 1361, Some(1371)), ("chr1", 1108, Some(1118)), ("chr2", 1246, None), ("chr2", 1239, Some(1249)), ("chr1", 1386, None), ("chr2", 1265, Some(1275)), ("chr2", 1479, None), ("chr1", 1669, Some(1679)), ("chr1", 1613, None), ("chr2", 1561, None), ("chr1", 1107, None), ("chr2", 1083, Some(1093)), ("chr1", 992, Some(1002)), ("chr1", 1178, Some(1188)), ("chr2", 1420, None), ("chr2", 1356, None), ("chr1", 1102, None), ("chr2", 1106, Some(1116)), ("chr1", 1225, None), ("chr2", 1343, None), ("chr2", 1290, None), ("chr2", 1107, Some(1117)), ("chr2", 1130, Some(1140)), ("chr1", 1336, Some(1346)), ("chr1", 1359, Some(1369)), ("chr1", 1165, None), ("chr1", 1680, None), ("chr2", 1405, Some(1415)), ("chr1", 1458, Some(1468)), ("chr2", 1049, Some(1059)), ("chr1", 1328, None), ("chr1", 1508, None), ("chr1", 1544, Some(1554)), ("chr1", 1180, Some(1190)), ("chr1", 1627, None), ("chr2", 1297, None), ("chr1", 1494, None), ("chr1", 1127, Some(1137)), ("chr2", 1684, Some(1694)), ("chr1", 1677, Some(1687)), ("chr2", 1570, Some(1580)), ("chr2", 1495, Some(1505)), ("chr1", 1060, None), ("chr1", 1082, None), ("chr1", 1606, None), ("chr1", 1273, None), ("chr1", 1181, None), ("chr1", 1016, None), ("chr1", 1440, Some(1450)), ("chr2", 1535, Some(1545)), ("chr1", 1197, Some(1207)), ("chr2", 1607, Some(1617)), ("chr2", 1371, Some(1381)), ("chr1", 1023, None), ("chr1", 1041, Some(1051)), ("chr2", 1677, Some(1687)), ("chr1", 1209, None), ("chr2", 1617, None), ("chr2", 1566, None), ("chr1", 1283, Some(1293)), ("chr1", 1677, None), ("chr1", 1391, None), ("chr2", 1385, None), ("chr2", 1201, None), ("chr1", 1277, Some(1287)), ("chr1", 1136, None)], 5, &[3, 12, 42, 45, 50, 52, 62, 69], 72, 10);
        check_rf(&[("chr1", 1219, Some(1229)), ("chr2", 1367, Some(1377)), ("chr2", 1503, Some(1513)), ("chr2", 1400, Some(1410)), ("chr2", 1508, None), ("chr2", 1160, None), ("chr1", 1435, None), ("chr2", 1508, Some(1518)), ("chr1", 1669, None), ("chr2", 1584, None), ("chr1", 1260, Some(1270)), ("chr1", 1483, Some(1493)), ("chr1", 1313, Some(1323)), ("chr1", 1299, Some(1309)), ("chr1", 1565, None), ("chr1", 1164, Some(1174)), ("chr1", 1684, Some(1694)), ("chr2", 1447, None), ("chr2", 1193, None), ("chr1", 1411, Some(1421)), ("chr1", 1315, Some(1325)), ("chr2", 1387, Some(1397)), ("chr1", 1187, None), ("chr1", 1155, Some(1165)), ("chr1", 1450, None), ("chr2", 1091, Some(1101)), ("chr1", 1444, Some(1454)), ("chr1", 1540, Some(1550)), ("chr2", 1645, None), ("chr1", 1576, Some(1586)), ("chr2", 1370, Some(1380)), ("chr1", 1648, Some(1658)), ("chr2", 1147, Some(1157)), ("chr1", 1061, None), ("chr1", 1192, Some(1202)), ("chr1", 1251, Some(1261)), ("chr1", 1259, Some(1269)), ("chr2", 990, Some(1000)), ("chr2", 1218, Some(1228)), ("chr2", 1116, None), ("chr1", 1107, Some(1117)), ("chr1", 1463, None), ("chr1", 1220, Some(1230)), ("chr1", 1320, None), ("chr2", 1667, None), ("chr2", 1219, Some(1229)), ("chr1", 1633, None), ("chr2", 1683, Some(1693)), ("chr2", 1271, Some(1281)), ("chr2", 1407, Some(1417)), ("chr1", 1563, Some(1573)), ("chr1", 1570, None), ("chr1", 1081, None), ("chr2", 999, Some(1009)), ("chr2", 1637, Some(1647)), ("chr1", 1481, None), ("chr2", 1444, None), ("chr1", 1146, None), ("chr1", 1663, Some(1673)), ("chr2", 1452, None), ("chr1", 1334, Some(1344)), ("chr1", 1676, Some(1686)), ("chr1", 1292, Some(1302)), ("chr1", 1108, Some(1118)), ("chr1", 1296, Some(1306)), ("chr2", 1183, None), ("chr2", 1115, Some(1125)), ("chr2", 1297, Some(1307)), ("chr2", 1382, Some(1392)), ("chr1", 1389, None), ("chr2", 1491, None), ("chr2", 1402, None), ("chr2", 1112, None), ("chr2", 1268, None), ("chr1", 1147, Some(1157)), ("chr2", 1360, Some(1370)), ("chr1", 1425, Some(1435)), ("chr2", 1233, Some(1243)), ("chr1", 1416, Some(1426)), ("chr2", 1679, None)], 4, &[37, 53], 78, 12);
        check_rf(&[("chr2", 1054, Some(1064)), ("chr1", 1404, Some(1414)), ("chr1", 1374, Some(1384)), ("chr2", 1080, Some(1090)), ("chr1", 1505, Some(1515)), ("chr1", 1092, None), ("chr1", 1130, Some(1140)), ("chr2", 1441, Some(1451)), ("chr1", 1320, None), ("chr1", 1576, None), ("chr1", 1556, None), ("chr2", 1312, None), ("chr1", 1092, Some(1102)), ("chr1", 1601, None), ("chr2", 1683, Some(1693)), ("chr1", 1578, Some(1588)), ("chr2", 1467, None), ("chr2", 1553, Some(1563)), ("chr1", 1096, None), ("chr2", 1338, Some(1348)), ("chr1", 1324, None), ("chr2", 1305, Some(1315)), ("chr2", 1515, Some(1525)), ("chr1", 1444, None), ("chr2", 1626, None), ("chr1", 1560, None), ("chr1", 1071, Some(1081)), ("chr1", 1359, Some(1369)), ("chr1", 1398, Some(1408)), ("chr1", 1307, None), ("chr1", 1178, Some(1188)), ("chr2", 1044, None), ("chr2", 1400, None), ("chr1", 1383, None), ("chr2", 1412, None), ("chr2", 1517, Some(1527)), ("chr2", 1144, None), ("chr2", 1559, Some(1569)), ("chr1", 1246, None), ("chr1", 1571, None), ("chr2", 1274, None), ("chr2", 1476, Some(1486)), ("chr2", 1188, None), ("chr1", 1164, Some(1174)), ("chr1", 1479, Some(1489)), ("chr1", 1527, Some(1537)), ("chr1", 1336, None), ("chr2", 1299, Some(1309)), ("chr1", 1404, None), ("chr2", 1215, Some(1225)), ("chr1", 1451, None), ("chr2", 990, Some(1000)), ("chr1", 1402, Some(1412)), ("chr1", 1597, Some(1607)), ("chr2", 1595, None), ("chr1", 1242, Some(1252)), ("chr1", 1049, Some(1059)), ("chr1", 1127, None), ("chr1", 1643, None), ("chr2", 1552, Some(1562)), ("chr2", 1158, None), ("chr2", 1262, Some(1272)), ("chr1", 1185, Some(1195)), ("chr2", 1603, Some(1613)), ("chr1", 1296, None), ("chr2", 1048, None), ("chr2", 1509, Some(1519)), ("chr2", 1260, None), ("chr2", 1668, Some(1678)), ("chr2", 1526, Some(1536)), ("chr2", 1478, None), ("chr2", 1112, None), ("chr1", 1106, Some(1116)), ("chr1", 1204, Some(1214)), ("chr2", 1673, Some(1683)), ("chr2", 1446, Some(1456)), ("chr1", 1646, None), ("chr1", 1168, Some(1178)), ("chr1", 1453, Some(1463)), ("chr2", 1421, None)], 2, &[], 80, 14);
        check_rf(&[("chr2", 1237, None), ("chr1", 1590, Some(1600)), ("chr2", 1312, None), ("chr2", 1394, Some(1404)), ("chr1", 1320, Some(1330)), ("chr2", 1197, Some(1207)), ("chr1", 1017, Some(1027)), ("chr1", 1433, None), ("chr2", 1504, None), ("chr1", 1589, Some(1599)), ("chr2", 1012, None), ("chr2", 1419, None), ("chr1", 1511, None), ("chr2", 1096, Some(1106)), ("chr2", 1664, Some(1674)), ("chr2", 1516, None), ("chr1", 1603, Some(1613)), ("chr2", 1046, None), ("chr2", 1680, None), ("chr1", 1485, None), ("chr1", 1359, None), ("chr2", 1095, Some(1105)), ("chr1", 1329, Some(1339)), ("chr1", 1570, Some(1580)), ("chr1", 1430, Some(1440)), ("chr1", 1638, None), ("chr2", 1649, None), ("chr1", 1296, None), ("chr1", 1141, Some(1151)), ("chr1", 1240, Some(1250)), ("chr1", 1258, Some(1268)), ("chr2", 1231, None), ("chr2", 1629, None), ("chr1", 1243, Some(1253)), ("chr2", 1344, Some(1354)), ("chr2", 1179, None), ("chr2", 1369, Some(1379)), ("chr2", 1494, Some(1504)), ("chr1", 1167, None), ("chr1", 1406, None), ("chr1", 1324, Some(1334)), ("chr1", 1378, Some(1388)), ("chr2", 1554, None), ("chr1", 1470, None), ("chr1", 1498, None), ("chr1", 1221, None), ("chr1", 1304, Some(1314)), ("chr1", 1055, None), ("chr2", 1099, None), ("chr1", 1448, None), ("chr2", 1132, None), ("chr2", 1689, Some(1699)), ("chr1", 1692, None), ("chr2", 1257, None), ("chr1", 1138, None), ("chr2", 1145, Some(1155)), ("chr1", 1689, None), ("chr1", 1313, None), ("chr2", 1320, Some(1330)), ("chr2", 1684, Some(1694)), ("chr1", 1443, Some(1453)), ("chr2", 1411, None), ("chr1", 1202, Some(1212)), ("chr1", 1093, None), ("chr2", 1075, Some(1085)), ("chr2", 1669, Some(1679)), ("chr2", 1630, Some(1640)), ("chr1", 1396, Some(1406)), ("chr1", 1601, None), ("chr2", 1464, None), ("chr1", 1384, Some(1394)), ("chr2", 1278, Some(1288)), ("chr1", 1209, Some(1219)), ("chr2", 1279, Some(1289)), ("chr2", 1300, Some(1310)), ("chr1", 1111, Some(1121)), ("chr2", 1427, Some(1437)), ("chr2", 1394, Some(1404)), ("chr1", 1513, Some(1523)), ("chr1", 1492, None)], 4, &[6], 79, 13);
        check_rf(&[("chr2", 1646, Some(1656)), ("chr1", 1401, None), ("chr2", 1420, None), ("chr2", 1322, Some(1332)), ("chr2", 1294, None), ("chr2", 1626, None), ("chr1", 1177, None), ("chr1", 1423, None), ("chr1", 1281, None), ("chr2", 1383, None), ("chr1", 1437, None), ("chr2", 1420, None), ("chr1", 1091, None), ("chr1", 1307, Some(1317)), ("chr2", 1380, None), ("chr2", 1638, Some(1648)), ("chr2", 1100, Some(1110)), ("chr2", 1521, Some(1531)), ("chr2", 1423, None), ("chr2", 1637, Some(1647)), ("chr2", 1327, Some(1337)), ("chr2", 1505, None), ("chr1", 1511, None), ("chr1", 1667, Some(1677)), ("chr1", 1047, Some(1057)), ("chr1", 1210, Some(1220)), ("chr2", 1452, None), ("chr2", 1535, Some(1545)), ("chr1", 1176, None), ("chr1", 1379, Some(1389)), ("chr2", 1360, Some(1370)), ("chr2", 1660, Some(1670)), ("chr1", 1034, Some(1044)), ("chr2", 1035, None), ("chr2", 1640, None), ("chr2", 1456, None), ("chr2", 1180, Some(1190)), ("chr1", 1464, None), ("chr2", 1137, None), ("chr2", 1565, Some(1575)), ("chr2", 1371, None), ("chr1", 1654, None), ("chr2", 1383, Some(1393)), ("chr2", 1003, Some(1013)), ("chr2", 1647, None), ("chr2", 1501, Some(1511)), ("chr1", 1295, Some(1305)), ("chr2", 1574, None), ("chr2", 1377, Some(1387)), ("chr1", 1589, None), ("chr1", 1603, None), ("chr2", 1645, None), ("chr2", 1499, Some(1509)), ("chr1", 1501, None), ("chr2", 1213, Some(1223)), ("chr1", 1255, Some(1265)), ("chr1", 1651, None), ("chr1", 1291, None), ("chr1", 1202, None), ("chr2", 1341, Some(1351)), ("chr2", 1356, None), ("chr1", 1581, Some(1591)), ("chr2", 1282, Some(1292)), ("chr1", 1088, Some(1098)), ("chr2", 1276, None), ("chr2", 1290, None), ("chr2", 1261, None), ("chr1", 1217, Some(1227)), ("chr2", 1203, None), ("chr2", 1340, Some(1350)), ("chr2", 1278, Some(1288)), ("chr2", 1130, Some(1140)), ("chr1", 1643, Some(1653)), ("chr1", 1510, Some(1520)), ("chr2", 1088, None), ("chr2", 1500, Some(1510)), ("chr1", 1341, Some(1351)), ("chr2", 1565, Some(1575)), ("chr2", 1337, None), ("chr1", 1264, None)], 2, &[], 80, 14);
    }
}
