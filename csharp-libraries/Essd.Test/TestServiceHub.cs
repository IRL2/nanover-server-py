using System;
using System.Collections.Generic;
using Newtonsoft.Json;
using NUnit.Framework;

namespace Essd.Test
{
    internal class TestServiceHub
    {

        [Test]
        public void TestMissingName()
        {
            try
            {
                var hub = new ServiceHub("{address:'1.2.3.4'}");
                
            }
            catch (ArgumentException)
            {
                Assert.Pass();
            }
            Assert.Fail();
        }
        
        [Test]
        public void TestMissingAddress()
        {
            try
            {
                var hub = new ServiceHub("{name:'test'}");
                
            }
            catch (ArgumentException)
            {
                Assert.Pass();
            }
            Assert.Fail();
        }
        
        [Test]
        public void TestValidJson()
        {
            var hub = new ServiceHub("{name:'test',address:'1.2.3.4'}");
            Assert.AreEqual("test", hub.Name);
            Assert.AreEqual("1.2.3.4", hub.Address);
        }
        
        [Test]
        public void TestInvalidJson()
        {
            try
            {
                var hub = new ServiceHub("{name:'test',address='1.2.3.4'");
            }
            catch (JsonReaderException)
            {
                Assert.Pass();
            }

            Assert.Fail();
        }

        [Test]
        public void TestToJson()
        {
            string expectedJson = "{'address':'1.2.3.4','name':'test','services':{'imd':54321},'id':'1234','essd_version':'1.0.0'}";
            
            var hub = new ServiceHub(expectedJson);
            var json = hub.ToJson();
            json = json.Replace('"', '\'');
            Assert.AreEqual(expectedJson, json);
        }

        [Test]
        public void TestEquality()
        {
            string json = "{'address':'1.2.3.4','name':'test','services':{'trajectory':54322,'imd':54321},'id':'1234','essd_version':'1.0.0'}";
            var hub = new ServiceHub(json);
            var secondHub = new ServiceHub("test", "1.2.3.4", "1234","1.0.0");
            var serviceDict = new Dictionary<string, int>();
            serviceDict["imd"] = 12345;
            serviceDict["trajectory"] = 54322;
            secondHub.Properties["services"] = serviceDict;
            
            Assert.AreEqual(hub, secondHub);
        }

        [Test]
        public void TestEqualityReference()
        {
            var hub = new ServiceHub("test", "1.2.3.4");
            var hub2 = hub; 
            
            Assert.IsTrue(hub.Equals(hub2));
        }
        
        [Test]
        public void TestEqualityNull()
        {
            ServiceHub hub = null;
            Assert.AreEqual(hub, null);
        }

        [Test]
        public void TestInEqualityNull()
        {
            var hub = new ServiceHub("test", "1.2.3.4");
            ServiceHub hub2 = null;
            Assert.IsFalse(hub.Equals(hub2));
        }
        
        [Test]
        public void TestObjectEqualityReference()
        {
            var hub = new ServiceHub("test", "1.2.3.4");
            var hub2 = (object) hub;
            Assert.IsTrue(hub.Equals(hub2));
            
        }

        [Test]
        public void TestObjectEqualityDifferentType()
        {
            ServiceHub hub = new ServiceHub("test", "1.2.3.4");
            int hub2 = 5;
            Assert.IsFalse(hub.Equals(hub2));
        }

        [Test]
        public void TestObjectInequalityNull()
        {
            var hub = new ServiceHub("test", "1.2.3.4");
            var hub2 = (object) null;
            Assert.IsFalse(hub.Equals(hub2));
        }
        
        [Test]
        public void TestObjectEquality()
        {
            string json = "{'address':'1.2.3.4','name':'test','id':'1234','services':{'trajectory':54322,'imd':54321}}";
            var hub = new ServiceHub(json);
            var secondHub = new ServiceHub("test", "1.2.3.4", "1234");
            var serviceDict = new Dictionary<string, int>();
            serviceDict["imd"] = 12345;
            serviceDict["trajectory"] = 54322;
            secondHub.Properties["services"] = serviceDict;

            var objectHub = (object) secondHub;
            Assert.IsTrue(hub.Equals(objectHub));
        }
        
        [Test]
        public void TestInequalityName()
        {
            string json = "{'address':'1.2.3.4','name':'test','id':'1234','services':{'trajectory':54322,'imd':54321}}";
            var hub = new ServiceHub(json);
            var secondHub = new ServiceHub("other", "1.2.3.4", id:"1234");
            var serviceDict = new Dictionary<string, int>();
            serviceDict["imd"] = 54322;
            serviceDict["trajectory"] = 54322;
            secondHub.Properties["services"] = serviceDict;
            
            Assert.AreNotEqual(hub, secondHub);
        }
        
        [Test]
        public void TestInequalityAddress()
        {
            string json = "{'address':'4.3.2.1','name':'test','services':{'trajectory':54322,'imd':54321}}";
            var hub = new ServiceHub(json);
            var secondHub = new ServiceHub("other", "1.2.3.4");
            var serviceDict = new Dictionary<string, int>();
            serviceDict["imd"] = 54322;
            serviceDict["trajectory"] = 54322;
            secondHub.Properties["services"] = serviceDict;
            
            Assert.AreNotEqual(hub, secondHub);
        }

        [Test]
        public void TestToString()
        {
            var hub = new ServiceHub("test", "1.2.3.4", id:"1234");
            Assert.AreEqual("test:1.2.3.4:1234", hub.ToString());

        }
        
        [Test]
        public void TestHashCode()
        {
            var hub = new ServiceHub("test", "1.2.3.4", id:"1234");
            int hash = "1234".GetHashCode();
            Assert.AreEqual(hash, hub.GetHashCode());
        }
    }
}