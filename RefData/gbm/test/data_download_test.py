import unittest, os, shutil
from gbm.GBMDataFinder import TriggerDataFinderFTP


class TestDataDownload(unittest.TestCase):
    
    def test_find_files(self):
        # test finding files on HEASARC
        download_dir = './downloads'
        dataFinder = TriggerDataFinderFTP(download_dir=download_dir)
        
        ctime = dataFinder.ctime_files(170313125, dets='n3')
        the_file = '/fermi/data/gbm/triggers/2017/bn170313125/current/glg_ctime_n3_bn170313125_v00.pha'
        self.assertEqual(ctime[0], the_file)
        
        cspec = dataFinder.cspec_files(170313125, dets=['n4'])
        the_file = '/fermi/data/gbm/triggers/2017/bn170313125/current/glg_cspec_n4_bn170313125_v00.pha'
        self.assertEqual(cspec[0], the_file)
    
        tte = dataFinder.tte_files(170313125, dets=['b0', 'n7'])
        the_file1 = '/fermi/data/gbm/triggers/2017/bn170313125/current/glg_tte_b0_bn170313125_v00.fit'
        the_file2 = '/fermi/data/gbm/triggers/2017/bn170313125/current/glg_tte_n7_bn170313125_v00.fit'
        self.assertEqual(tte[0], the_file1)
        self.assertEqual(tte[1], the_file2)


        ctime_rsp1 = dataFinder.rsp_files(170313125, ctime=True, dets='n3')
        the_file = '/fermi/data/gbm/triggers/2017/bn170313125/current/glg_ctime_n3_bn170313125_v01.rsp'
        self.assertEqual(ctime_rsp1[0], the_file)

        cspec_rsp1 = dataFinder.rsp_files(170313125, cspec=True, dets=['n4'])
        the_file = '/fermi/data/gbm/triggers/2017/bn170313125/current/glg_cspec_n4_bn170313125_v01.rsp'
        self.assertEqual(cspec_rsp1[0], the_file)
        
        rsp2s = dataFinder.rsp_files(170313125, cspec=True, rsp2=True, dets='b0')
        the_file1 = '/fermi/data/gbm/triggers/2017/bn170313125/current/glg_cspec_b0_bn170313125_v01.rsp2'
        self.assertEqual(rsp2s[0], the_file1)
        
        trigdat = dataFinder.trigdat_file(170313125)
        the_file = '/fermi/data/gbm/triggers/2017/bn170313125/current/glg_trigdat_all_bn170313125_v01.fit'
        self.assertEqual(trigdat, the_file)
        
    def test_not_exist(self):
        # test non-existent trigger number
        download_dir = './downloads'
        dataFinder = TriggerDataFinderFTP(download_dir=download_dir)
        
        ctime = dataFinder.ctime_files(13141516)
        self.assertEqual(len(ctime), 0)
    
    def test_download_files(self):
        # test the download of files
        download_dir = './downloads'
        dataFinder = TriggerDataFinderFTP(download_dir=download_dir)
        
        ctime = dataFinder.ctime_files(170313125, dets='n0')
        the_file = dataFinder.ftp_get(ctime)
        self.assertTrue(os.path.isfile(the_file[0]))
        try:
            shutil.rmtree(download_dir)
        except:
            pass

if __name__ == '__main__':
    unittest.main()
